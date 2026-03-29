#!/usr/bin/env python3
"""
SI §4.3 — κ fine-tuning on compact2 with patristic distance regression.

Loads compact2 (V15Model), freezes ALL parameters except the curvature
scalar k = encoder.manifold.k, then runs grouped batches and minimises
MSE between Poincaré distances and taxonomic rank distances.

Key design decisions vs the BiosphereCodec attempt:
  - Embeddings are computed via encode_angular_only() and DETACHED from k.
    Gradient flows only through the distance formula d(z_i, z_j; k).
    This is identical to the telescope sweep (fit_kappa_telescope.py) but
    done online in a training loop rather than post-hoc.
  - Rank distances are batch-normalised to match the current mean distance
    scale, so the loss is purely about relative ordering (same-genus closer
    than cross-domain), not absolute magnitude.
  - No warmup needed: compact2 embeddings are already structured.
  - GPU safety: checks free memory before each batch, yields to other
    processes if >80% VRAM is used.

SI prediction: κ converges to 1.247 ± 0.003.
Current compact2 κ: 1.2505.

Usage (on inference server):
    python finetune_kappa_patristic.py \\
        --checkpoint /fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt \\
        --manifest /fast/sentrybio/data/manifest_local.csv \\
        --output-dir ./kappa_finetune \\
        --steps 3000

    # 5-seed sweep
    python finetune_kappa_patristic.py ... --sweep
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import random
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Sampler


# ── V15Model import ───────────────────────────────────────────────────────
sys.path.insert(0, "/fast/sentrybio/scripts")
sys.path.insert(0, str(Path.home()))
from train_v15_5_phase2c_next import V15Model


# ── Poincaré distance (differentiable w.r.t. c) ──────────────────────────

def poincare_dist_mat(Z: torch.Tensor, c: torch.Tensor, eps: float = 1e-7) -> torch.Tensor:
    """Full pairwise Poincaré distance matrix. Gradient flows through c."""
    B = Z.shape[0]
    u = Z.unsqueeze(1).expand(B, B, -1)
    v = Z.unsqueeze(0).expand(B, B, -1)
    diff_sq = ((u - v) ** 2).sum(-1)
    u_sq    = (u ** 2).sum(-1)
    v_sq    = (v ** 2).sum(-1)
    denom   = ((1 - c * u_sq) * (1 - c * v_sq)).clamp(min=eps)
    arg     = (1 + 2 * c * diff_sq / denom).clamp(min=1 + eps)
    return torch.acosh(arg) / torch.sqrt(c + eps)


# ── Taxonomy rank distances ───────────────────────────────────────────────

_RANK_COLS = ["genus", "family", "order", "class", "phylum", "domain"]
# Distances on a 7-point scale: same genus=1/7, cross-domain=7/7=1.0
# (never zero — avoids distance-collapse gradient)
_RANK_VALS = [1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1.0]


def tax_rank_dist(row_i: dict, row_j: dict) -> float:
    for level, col in enumerate(_RANK_COLS):
        vi = row_i.get(col, "").strip()
        vj = row_j.get(col, "").strip()
        if vi and vj and vi == vj:
            return _RANK_VALS[level]
    return _RANK_VALS[6]


def build_rank_matrix(tax_rows: List[dict], device: torch.device) -> torch.Tensor:
    B = len(tax_rows)
    D = torch.zeros(B, B, dtype=torch.float32)
    for i in range(B):
        for j in range(i + 1, B):
            d = tax_rank_dist(tax_rows[i], tax_rows[j])
            D[i, j] = D[j, i] = d
    return D.to(device)


# ── Dataset ───────────────────────────────────────────────────────────────

class GenomeDataset(Dataset):
    def __init__(
        self,
        manifest_path: str,
        max_len: int = 8192,
        max_genomes: Optional[int] = None,
        min_genus_count: int = 4,
        vocab_size: int = 4096,
    ):
        rows = []
        with open(manifest_path, newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                tp = (row.get("tokenized_path") or "").strip()
                g  = (row.get("genus") or "").strip()
                if tp and g and os.path.exists(tp):
                    rows.append(row)

        from collections import Counter
        gc = Counter(r["genus"] for r in rows)
        valid = {g for g, cnt in gc.items() if cnt >= min_genus_count}
        rows = [r for r in rows if r["genus"].strip() in valid]

        if max_genomes and len(rows) > max_genomes:
            random.Random(42).shuffle(rows)
            rows = rows[:max_genomes]

        self.rows = rows
        genera = sorted(set(r["genus"].strip() for r in rows))
        self.genus_to_id = {g: i for i, g in enumerate(genera)}
        self.labels = [self.genus_to_id[r["genus"].strip()] for r in rows]
        self.max_len = max_len
        self.vocab_size = vocab_size
        print(f"Dataset: {len(rows)} genomes, {len(genera)} genera (≥{min_genus_count}/genus)")

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, idx):
        row = self.rows[idx]
        tokens = np.load(row["tokenized_path"]).astype(np.int64)
        tokens = np.clip(tokens, 0, self.vocab_size - 1)
        if len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start : start + self.max_len]
        return torch.from_numpy(tokens), self.labels[idx], idx


# ── Grouped sampler ───────────────────────────────────────────────────────

class GenusGroupedBatchSampler(Sampler[List[int]]):
    def __init__(self, labels, batch_size=16, members_per_genus=4, seed=42):
        self.members_per_genus = members_per_genus
        self.genera_per_batch  = batch_size // members_per_genus
        self.rng = np.random.RandomState(seed)

        gidx: Dict[int, List[int]] = defaultdict(list)
        for i, lab in enumerate(labels):
            gidx[lab].append(i)
        self.gidx = {g: v for g, v in gidx.items() if len(v) >= members_per_genus}
        self.gids = list(self.gidx.keys())
        total = sum(len(v) for v in self.gidx.values())
        self._len = max(1, total // batch_size)
        print(f"  Sampler: {self.genera_per_batch} genera × {members_per_genus} "
              f"= {batch_size}/batch, {len(self.gids)} eligible genera")

    def __iter__(self):
        order = self.gids.copy()
        self.rng.shuffle(order)
        pools = {g: self.gidx[g].copy() for g in order}
        for g in pools: self.rng.shuffle(pools[g])
        g_ptr = 0
        for _ in range(self._len):
            batch = []
            for _ in range(self.genera_per_batch):
                g = order[g_ptr % len(order)]
                g_ptr += 1
                pool = pools[g]
                for _ in range(self.members_per_genus):
                    if not pool:
                        pool = self.gidx[g].copy()
                        self.rng.shuffle(pool)
                        pools[g] = pool
                    batch.append(pool.pop())
            yield batch

    def __len__(self):
        return self._len


def collate_fn(batch):
    tokens_list, labels, idxs = zip(*batch)
    max_len = max(t.shape[0] for t in tokens_list)
    padded = torch.zeros(len(batch), max_len, dtype=torch.long)
    for i, t in enumerate(tokens_list):
        padded[i, :t.shape[0]] = t
    return padded, torch.tensor(labels, dtype=torch.long), list(idxs)


# ── Load compact2 ─────────────────────────────────────────────────────────

def load_compact2(checkpoint_path: str, device: torch.device) -> V15Model:
    ckpt  = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = ckpt.get("model_state_dict", ckpt)

    latent_dim = state["ode_flow.field.0.weight_orig"].shape[1] - 1
    ode_hidden  = state["ode_flow.field.0.bias"].shape[0]
    counts = {
        "bact_fam": state["bact_fam.prototypes"].shape[0],
        "arch_fam": state["arch_fam.prototypes"].shape[0],
        "euk_fam":  state["euk_fam.prototypes"].shape[0],
        "bact_gen": state["bact_gen.prototypes"].shape[0],
        "arch_gen": state["arch_gen.prototypes"].shape[0],
        "euk_gen":  state["euk_gen.prototypes"].shape[0],
    }
    vocab_size = state["encoder.encoder.embed.weight"].shape[0]

    model = V15Model(
        vocab_size=vocab_size, latent_dim=latent_dim,
        counts=counts, ode_hidden=ode_hidden,
    ).to(device)

    state_filtered = {k: v for k, v in state.items() if "curvature_history" not in k}
    missing, unexpected = model.load_state_dict(state_filtered, strict=False)
    if missing:    print(f"  Missing keys:    {len(missing)}")
    if unexpected: print(f"  Unexpected keys: {len(unexpected)}")

    print(f"  compact2 loaded. live_kappa = {model.live_kappa:.6f}")
    return model


# ── GPU memory guard ──────────────────────────────────────────────────────

def gpu_free_fraction() -> float:
    """Returns fraction of GPU memory that is free (0–1). Returns 1.0 on CPU."""
    if not torch.cuda.is_available():
        return 1.0
    free, total = torch.cuda.mem_get_info()
    return free / total


# ── Training ──────────────────────────────────────────────────────────────

def train(args) -> float:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[
            logging.FileHandler(outdir / "training.log"),
            logging.StreamHandler(),
        ],
    )
    log = logging.getLogger("kappa_finetune")

    # ── Check GPU before loading model ────────────────────────────────────
    if torch.cuda.is_available():
        free_frac = gpu_free_fraction()
        log.info(f"GPU: {free_frac*100:.1f}% free before model load")
        if free_frac < 0.25:
            log.warning("GPU <25% free — falling back to CPU to avoid disrupting other processes")
            device = torch.device("cpu")

    # ── Load compact2 ─────────────────────────────────────────────────────
    log.info(f"Loading compact2 from {args.checkpoint}")
    model = load_compact2(args.checkpoint, device)
    kappa_init = model.live_kappa

    # ── Freeze everything except curvature scalar k ───────────────────────
    k_param = model.encoder.manifold.k
    for name, param in model.named_parameters():
        param.requires_grad_(False)
    k_param.requires_grad_(True)
    n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    log.info(f"Frozen all params. Trainable: {n_trainable} (curvature scalar only)")
    log.info(f"Initial κ = {kappa_init:.6f}  (compact2 training κ)")
    log.info(f"SI §4.3 prediction: κ → 1.247 ± 0.003")

    model.eval()  # inference mode — no dropout, no batch-norm updates

    # ── Dataset ───────────────────────────────────────────────────────────
    dataset = GenomeDataset(
        args.manifest,
        max_len=args.max_len,
        max_genomes=args.max_genomes,
        min_genus_count=args.members_per_genus,
    )
    sampler = GenusGroupedBatchSampler(
        dataset.labels,
        batch_size=args.batch_size,
        members_per_genus=args.members_per_genus,
        seed=args.seed,
    )
    loader = DataLoader(dataset, batch_sampler=sampler, collate_fn=collate_fn, num_workers=0)

    # ── Optimizer: only k_param ───────────────────────────────────────────
    # Very small lr: compact2 embeddings are calibrated for κ≈1.25.
    # We only want to probe whether patristic pushes κ slightly up or down
    # from that equilibrium — not let it drift far enough to enter a cliff regime.
    optimizer = torch.optim.Adam([k_param], lr=args.lr)

    kappa_history = []
    loader_iter   = iter(loader)
    step          = 0

    while step < args.steps:
        # GPU safety: pause if another process is using >80% VRAM
        if torch.cuda.is_available():
            free_frac = gpu_free_fraction()
            if free_frac < 0.20:
                log.info(f"Step {step}: GPU <20% free ({free_frac*100:.1f}%), waiting 30s...")
                torch.cuda.empty_cache()
                time.sleep(30)
                continue

        try:
            tokens, tax_ids, batch_idxs = next(loader_iter)
        except StopIteration:
            loader_iter = iter(loader)
            tokens, tax_ids, batch_idxs = next(loader_iter)

        tokens = tokens.to(device)

        # ── Extract embeddings (no gradient — frozen encoder) ─────────────
        with torch.no_grad():
            z = model.encode_angular_only(tokens)  # [B, emb_dim]

        # ── Build [B,B] rank distance matrix ─────────────────────────────
        tax_rows = [dataset.rows[i] for i in batch_idxs]
        D_rank   = build_rank_matrix(tax_rows, device)  # values in [1/7, 1.0]

        # ── Poincaré distances at current κ (gradient flows through k) ────
        c   = k_param.abs() + 1e-4   # ensure positive curvature
        hyp = poincare_dist_mat(z.float(), c)

        # ── Batch-normalised patristic loss ───────────────────────────────
        # Scale rank targets to match mean actual distance — preserves
        # relative ordering without being sensitive to absolute scale.
        # V15 embeddings are already structured, so same-genus pairs are
        # genuinely closer than cross-domain pairs in each batch.
        mask = torch.triu(torch.ones_like(hyp, dtype=torch.bool), diagonal=1)
        mean_hyp  = hyp[mask].mean().detach()
        mean_rank = D_rank[mask].mean()
        scaled    = D_rank * (mean_hyp / mean_rank.clamp(min=1e-6))
        dist_loss = F.mse_loss(hyp[mask], scaled[mask])

        optimizer.zero_grad()
        dist_loss.backward()
        torch.nn.utils.clip_grad_norm_([k_param], 0.5)
        optimizer.step()

        # Clamp near compact2's training value — embeddings are only meaningful
        # in the geometry they were trained for (κ≈1.25). Straying far causes
        # ball-boundary cliffs and uninformative gradients.
        with torch.no_grad():
            k_param.clamp_(1.0, 2.0)

        kappa = model.live_kappa
        grad  = k_param.grad.item() if k_param.grad is not None else 0.0
        step += 1

        kappa_history.append({
            "step": step, "kappa": kappa,
            "dist_loss": dist_loss.item(),
            "mean_hyp": mean_hyp.item(),
            "c_grad": grad,
        })

        if step % args.log_every == 0:
            log.info(
                f"Step {step:>4d} | κ={kappa:.6f} | ∇κ={grad:+.4e} | "
                f"dist={dist_loss.item():.5f} | d̄={mean_hyp.item():.4f}"
            )

        if step % 500 == 0:
            torch.save({"step": step, "kappa": kappa, "k_param": k_param.item()},
                       outdir / f"kappa_{step}.pt")

    # ── Final result ──────────────────────────────────────────────────────
    final_kappa = model.live_kappa
    plateau = [h["kappa"] for h in kappa_history[len(kappa_history)*4//5:]]
    plateau_mean = float(np.mean(plateau))
    plateau_std  = float(np.std(plateau))

    theory = (1.6 * math.log(2)) ** 2
    log.info("=" * 60)
    log.info(f"COMPLETE (seed={args.seed})")
    log.info(f"  Initial κ:      {kappa_init:.6f}")
    log.info(f"  Final κ:        {final_kappa:.6f}")
    log.info(f"  Plateau κ:      {plateau_mean:.6f} ± {plateau_std:.6f}")
    log.info(f"  SI prediction:  1.247000 ± 0.003000")
    log.info(f"  Theory (h·ln2)²:{theory:.6f}")
    log.info(f"  Δ from SI pred: {abs(plateau_mean - 1.247) / 1.247 * 100:.1f}%")

    result = {
        "seed": args.seed,
        "kappa_init": kappa_init,
        "kappa_final": final_kappa,
        "kappa_plateau_mean": plateau_mean,
        "kappa_plateau_std": plateau_std,
        "theory_kappa": theory,
        "si_prediction": 1.247,
    }
    with open(outdir / "result.json", "w") as f:
        json.dump(result, f, indent=2)
    with open(outdir / "kappa_history.json", "w") as f:
        json.dump(kappa_history, f, indent=2)

    return plateau_mean


def main():
    p = argparse.ArgumentParser(description="SI §4.3: κ fine-tuning on compact2 with patristic")
    p.add_argument("--checkpoint", default="/fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt")
    p.add_argument("--manifest",   default="/fast/sentrybio/data/manifest_local.csv")
    p.add_argument("--output-dir", default="./kappa_finetune")
    p.add_argument("--max-len",    type=int, default=8192)
    p.add_argument("--max-genomes",type=int, default=None)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--members-per-genus", type=int, default=4)
    p.add_argument("--lr",         type=float, default=1e-4,
                   help="Learning rate for κ. Keep small (1e-4 to 1e-5) — compact2 "
                        "embeddings are calibrated for κ≈1.25, large steps cause cliffs.")
    p.add_argument("--steps",      type=int, default=3000)
    p.add_argument("--seed",       type=int, default=42)
    p.add_argument("--log-every",  type=int, default=10)
    p.add_argument("--sweep", action="store_true",
                   help="5-seed sweep (seeds 0,42,123,456,789) to reproduce SI Table 1")
    args = p.parse_args()

    if args.sweep:
        seeds   = [0, 42, 123, 456, 789]
        results = []
        base    = Path(args.output_dir)
        for seed in seeds:
            args.seed       = seed
            args.output_dir = str(base / f"seed_{seed}")
            kappa = train(args)
            results.append({"seed": seed, "kappa": kappa})
            with open(base / "sweep_results.json", "w") as f:
                json.dump(results, f, indent=2)

        kappas = [r["kappa"] for r in results]
        mean_k, std_k = float(np.mean(kappas)), float(np.std(kappas))
        theory = (1.6 * math.log(2)) ** 2
        print("\n" + "=" * 60)
        print("5-SEED SWEEP — SI §4.3")
        print("=" * 60)
        for r in results:
            print(f"  seed={r['seed']:>3d}  κ={r['kappa']:.6f}")
        print(f"\nκ = {mean_k:.6f} ± {std_k:.6f}  (CV={std_k/mean_k*100:.1f}%)")
        print(f"SI prediction: 1.247000 ± 0.003000  (CV=0.2%)")
        print(f"Theory:        {theory:.6f}")
        with open(base / "sweep_results.json", "w") as f:
            json.dump({
                "seeds": results,
                "kappa_mean": mean_k, "kappa_std": std_k,
                "kappa_cv_pct": std_k/mean_k*100,
                "theory_kappa": theory, "si_prediction": 1.247,
            }, f, indent=2)
    else:
        train(args)


if __name__ == "__main__":
    main()
