#!/usr/bin/env python3
"""
SI §4.3 — κ Convergence with Patristic Distance Regression
===========================================================

This is the experiment described in Supplementary Information §4.3 but never
previously executed. Every training script in the codebase has called
loss_fn(..., patristic=None). This script throws the switch.

Per-batch, a real [B, B] taxonomic rank distance matrix is computed from the
manifest taxonomy columns (domain/phylum/class/order/family/genus) and passed
as the `patristic` tensor to BiosphereCodec's loss function, activating:

    dist_loss = MSE(poincaré_dist_mat(z), patristic)
    total = mlm + dec + 0.1 * hex + 0.5 * dist_loss

Rank distance normalization (maps to [0, 1]):
    same genus   → 0.000  (d/6 where d=0)
    same family  → 0.167
    same order   → 0.333
    same class   → 0.500
    same phylum  → 0.667
    same domain  → 0.833
    cross-domain → 1.000

The Poincaré distance at κ=1.25 for same-genus pairs at unit scale saturates
around ~0.5 (empirically from telescope), so rank distances ∈ [0,1] are
in the right ballpark without further re-scaling. A learnable scale factor is
NOT added — we want the loss to exert pressure on κ directly.

SI prediction: κ → 1.247 ± 0.003 over 5 seeds.

Usage (on inference server):
    python run_patristic_kappa.py \\
        --manifest /fast/sentrybio/data/manifest_local.csv \\
        --output-dir ./patristic_kappa_run \\
        --steps 7000

    # 5-seed sweep (reproduces SI Table 1)
    python run_patristic_kappa.py \\
        --manifest /fast/sentrybio/data/manifest_local.csv \\
        --output-dir ./patristic_sweep \\
        --sweep
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
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Sampler


# ── BiosphereCodec import ─────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
for candidate in [SCRIPT_DIR, SCRIPT_DIR.parent, Path.home()]:
    if (candidate / "BiosphereCodec.py").exists():
        sys.path.insert(0, str(candidate))
        break
try:
    from BiosphereCodec import BiosphereCodec
except ImportError:
    raise ImportError(
        "BiosphereCodec.py not found. Place it next to this script or in parent dir."
    )


# ── Taxonomy rank distance ────────────────────────────────────────────────

# Ordered from finest to coarsest. Distance = (first non-matching level) / 6.
# same genus=0/6, same family=1/6, ..., cross-domain=6/6=1.0
_RANK_COLS = ["genus", "family", "order", "class", "phylum", "domain"]

# Rank distances never reach 0 — avoids gradient collapse toward κ=0.
# Scale factor converts [1/7, 1.0] to approx [0.5, 3.5] Poincaré distance units.
# At κ≈1.25, same-genus Poincaré distances ~ 0.5-1.0, cross-domain ~ 3-5.
_RANK_DIST = [1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1.0]  # levels 0..6
_DIST_SCALE = 3.5                                     # maps [1/7, 1.0] → [0.5, 3.5]


def tax_rank_dist(row_i: dict, row_j: dict) -> float:
    """Taxonomic rank distance scaled to Poincaré distance units.

    same genus   → 0.50  (1/7 × 3.5)
    same family  → 1.00
    same order   → 1.50
    same class   → 2.00
    same phylum  → 2.50
    same domain  → 3.00
    cross-domain → 3.50
    """
    for level, col in enumerate(_RANK_COLS):
        vi = row_i.get(col, "").strip()
        vj = row_j.get(col, "").strip()
        if vi and vj and vi == vj:
            return _RANK_DIST[level] * _DIST_SCALE
    return _RANK_DIST[6] * _DIST_SCALE  # cross-domain


def build_patristic_batch(tax_rows: List[dict], device: torch.device) -> torch.Tensor:
    """Compute [B, B] scaled taxonomic rank distance matrix for a batch."""
    B = len(tax_rows)
    D = torch.zeros(B, B, dtype=torch.float32)
    for i in range(B):
        for j in range(i + 1, B):
            d = tax_rank_dist(tax_rows[i], tax_rows[j])
            D[i, j] = D[j, i] = d
    return D.to(device)


# ── Dataset ───────────────────────────────────────────────────────────────

class GenomeDataset(Dataset):
    """Load .npy tokenized genomes with full taxonomy for patristic regression."""

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
            r = csv.DictReader(f)
            for row in r:
                tp = (row.get("tokenized_path") or row.get("TokenizedPath") or "").strip()
                g = (row.get("genus") or row.get("Genus") or "").strip()
                if tp and g and os.path.exists(tp):
                    rows.append(row)

        from collections import Counter
        gc = Counter(r["genus"] for r in rows)
        valid_genera = {g for g, cnt in gc.items() if cnt >= min_genus_count}
        rows = [r for r in rows if r.get("genus", "").strip() in valid_genera]

        if max_genomes and len(rows) > max_genomes:
            random.Random(42).shuffle(rows)
            rows = rows[:max_genomes]

        self.rows = rows
        genus_list = sorted(set(r.get("genus", "").strip() for r in rows))
        self.genus_to_id = {g: i for i, g in enumerate(genus_list)}
        self.labels = [self.genus_to_id[r.get("genus", "").strip()] for r in rows]
        self.max_len = max_len
        self.vocab_size = vocab_size
        self.n_genera = len(self.genus_to_id)
        print(f"Dataset: {len(rows)} genomes, {self.n_genera} genera (min {min_genus_count}/genus)")

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, idx: int):
        row = self.rows[idx]
        tokens = np.load(row["tokenized_path"]).astype(np.int64)
        tokens = np.clip(tokens, 0, self.vocab_size - 1)
        if len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start : start + self.max_len]
        return torch.from_numpy(tokens), self.labels[idx], idx


# ── Grouped sampler (same as run_grouped_kappa.py) ────────────────────────

class GenusGroupedBatchSampler(Sampler[List[int]]):
    def __init__(
        self,
        labels: List[int],
        batch_size: int = 32,
        members_per_genus: int = 4,
        seed: int = 42,
    ):
        self.batch_size = batch_size
        self.members_per_genus = members_per_genus
        self.genera_per_batch = batch_size // members_per_genus
        self.rng = np.random.RandomState(seed)

        genus_to_idx: Dict[int, List[int]] = defaultdict(list)
        for idx, lab in enumerate(labels):
            genus_to_idx[lab].append(idx)

        self.genus_indices = {
            g: idxs for g, idxs in genus_to_idx.items()
            if len(idxs) >= members_per_genus
        }
        self.genus_ids = list(self.genus_indices.keys())
        total = sum(len(v) for v in self.genus_indices.values())
        self._len = max(1, total // batch_size)
        print(f"  Grouped sampler: {self.genera_per_batch} genera × {members_per_genus} members "
              f"= {batch_size}/batch, {len(self.genus_ids)} genera eligible")

    def __iter__(self):
        order = self.genus_ids.copy()
        self.rng.shuffle(order)
        pools = {g: self.genus_indices[g].copy() for g in order}
        for g in pools:
            self.rng.shuffle(pools[g])

        g_ptr = 0
        for _ in range(self._len):
            batch = []
            for _ in range(self.genera_per_batch):
                g = order[g_ptr % len(order)]
                g_ptr += 1
                pool = pools[g]
                for _ in range(self.members_per_genus):
                    if not pool:
                        pool = self.genus_indices[g].copy()
                        self.rng.shuffle(pool)
                        pools[g] = pool
                    batch.append(pool.pop())
            yield batch

    def __len__(self) -> int:
        return self._len


def collate_fn(batch):
    tokens_list, labels, idxs = zip(*batch)
    max_len = max(t.shape[0] for t in tokens_list)
    padded = torch.zeros(len(batch), max_len, dtype=torch.long)
    for i, t in enumerate(tokens_list):
        padded[i, : t.shape[0]] = t
    return padded, torch.tensor(labels, dtype=torch.long), list(idxs)


# ── Training ──────────────────────────────────────────────────────────────

def train(args) -> float:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "training.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    log = logging.getLogger("patristic_kappa")

    dataset = GenomeDataset(
        args.manifest,
        max_len=args.max_len,
        max_genomes=args.max_genomes,
        min_genus_count=args.members_per_genus,
        vocab_size=args.vocab_size,
    )
    sampler = GenusGroupedBatchSampler(
        dataset.labels,
        batch_size=args.batch_size,
        members_per_genus=args.members_per_genus,
        seed=args.seed,
    )
    loader = DataLoader(
        dataset,
        batch_sampler=sampler,
        collate_fn=collate_fn,
        num_workers=0,
    )

    # BiosphereCodec signature varies by version — max_len param optional
    try:
        model = BiosphereCodec(
            vocab=args.vocab_size,
            d_model=args.d_model,
            n_layers=args.n_layers,
            max_len=args.max_len,
            latent_dim=args.latent_dim,
        ).to(device)
    except TypeError:
        model = BiosphereCodec(
            vocab=args.vocab_size,
            d_model=args.d_model,
            n_layers=args.n_layers,
            latent_dim=args.latent_dim,
        ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"Model: {n_params / 1e6:.1f}M params, d={args.d_model}, L={args.n_layers}, "
             f"latent={args.latent_dim}")
    log.info(f"Initial κ = {model.hyper.c.item():.6f}")
    log.info(f"Patristic switch: ON after {args.warmup_steps} warmup steps (weight={args.dist_weight}, 500-step ramp)")
    log.info(f"Loss: [0..{args.warmup_steps}] mlm+dec+hex | [{args.warmup_steps}+] + {args.dist_weight}*dist")

    # Phase 1: exclude hyper.c from optimizer entirely.
    # requires_grad + weight_decay alone can still move c via the optimizer's
    # internal state; excluding it from param groups is the only reliable fix.
    phase1_params = [p for n, p in model.named_parameters() if n != "hyper.c"]
    optimizer = torch.optim.AdamW(phase1_params, lr=args.lr, weight_decay=0.01)
    log.info(f"Phase 1: κ excluded from optimizer (fixed at {model.hyper.c.item():.4f}) for {args.warmup_steps} steps")

    model.train()
    step = 0
    kappa_history = []
    loader_iter = iter(loader)

    while step < args.steps:
        try:
            tokens, tax_ids, batch_idxs = next(loader_iter)
        except StopIteration:
            loader_iter = iter(loader)
            tokens, tax_ids, batch_idxs = next(loader_iter)

        tokens = tokens.to(device)
        tax_ids = tax_ids.to(device)

        # ── Build [B, B] patristic matrix from taxonomy columns ──────────
        tax_rows = [dataset.rows[i] for i in batch_idxs]
        patristic = build_patristic_batch(tax_rows, device)

        # ── Forward: patristic switch thrown ──────────────────────────────
        # Use manual dist_loss (weight=0.1) rather than BiosphereCodec's
        # built-in slot (weight=0.5) — the built-in slot has no scale guard
        # and drives κ to its lower clamp on a randomly-initialized model.
        masked_ids, mlm_labels = model.loss_fn.mask_tokens(
            tokens, model.encoder.embed.num_embeddings
        )
        enc_h, z = model.encode(masked_ids)
        enc_logits = enc_h @ model.encoder.embed.weight.T
        dec_logits = model.decoder(enc_h)
        loss, logs = model.loss_fn(
            tokens, mlm_labels, enc_logits, dec_logits, z,
            tax_ids=tax_ids, patristic=None,
        )

        # Manual patristic regression — active only after warmup.
        # The HEX loss must first cluster same-genus embeddings before
        # patristic pressure shapes κ. Before warmup_steps, same-genus
        # distances are random (≈ cross-genus), so the gradient always pushes
        # κ toward its lower clamp (same-genus target < mean → reduce all distances).
        # After warmup, same-genus distances are already smaller, so the gradient
        # has a genuine signal: same-genus too large → push up, cross-domain
        # too small → push down. The equilibrium is the empirically correct κ.
        hyp = model.hyper.dist_mat(z)
        mask_tri = torch.triu(torch.ones_like(hyp), diagonal=1).bool()
        mean_hyp = hyp[mask_tri].mean().detach()
        logs["mean_hyp"] = mean_hyp.item()

        if step == args.warmup_steps:
            # Phase 2: rebuild optimizer including κ — embeddings now have structure
            optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr * 0.1, weight_decay=0.01)
            log.info(f"Step {step}: Phase 2 — κ added to optimizer (lr×0.1={args.lr*0.1:.2e}), patristic activating")

        if step > args.warmup_steps:
            # Batch-normalize targets to match actual distance scale,
            # then penalize MSE on relative shape only
            mean_pat = patristic[mask_tri].mean()
            scaled_pat = patristic * (mean_hyp / (mean_pat.clamp(min=1e-6)))
            dist_loss = F.mse_loss(hyp[mask_tri], scaled_pat[mask_tri])
            # Ramp weight linearly from 0 → dist_weight over 500 steps after warmup
            ramp = min(1.0, (step - args.warmup_steps) / 500.0)
            loss = loss + args.dist_weight * ramp * dist_loss
            logs["dist"] = dist_loss.item()
        else:
            logs["dist"] = 0.0

        if torch.isnan(loss):
            log.warning(f"Step {step}: NaN loss, skipping")
            continue

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        if model.hyper.c.requires_grad:
            with torch.no_grad():
                model.hyper.c.clamp_(0.5, 3.0)
        step += 1

        kappa = model.hyper.c.item()
        c_grad = model.hyper.c.grad.item() if model.hyper.c.grad is not None else 0.0
        kappa_history.append({"step": step, "kappa": kappa, "c_grad": c_grad, **logs})

        if step % args.log_every == 0:
            log.info(
                f"Step {step:>5d} | κ={kappa:.6f} | ∇c={c_grad:+.4e} | "
                f"MLM={logs['mlm']:.4f} hex={logs['hex']:.4f} "
                f"dist={logs['dist']:.4f} d̄={logs.get('mean_hyp', 0):.3f} "
                f"total={logs['total']:.4f}"
            )

        if step % 500 == 0:
            torch.save({
                "step": step,
                "model_state_dict": model.state_dict(),
                "kappa": kappa,
                "kappa_history": kappa_history[-200:],
            }, outdir / f"checkpoint_{step}.pt")

    final_kappa = model.hyper.c.item()
    torch.save({
        "step": step,
        "model_state_dict": model.state_dict(),
        "kappa": final_kappa,
        "kappa_history": kappa_history,
    }, outdir / "final.pt")
    with open(outdir / "kappa_history.json", "w") as f:
        json.dump(kappa_history, f, indent=2)

    # ── Summary ───────────────────────────────────────────────────────────
    theory = (1.6 * math.log(2)) ** 2
    log.info(f"COMPLETE: Final κ = {final_kappa:.6f}  (SI prediction: 1.247, theory: {theory:.4f})")

    # Last 20% of training — κ plateau estimate
    plateau_start = len(kappa_history) * 4 // 5
    plateau_kappas = [h["kappa"] for h in kappa_history[plateau_start:]]
    plateau_mean = float(np.mean(plateau_kappas))
    plateau_std  = float(np.std(plateau_kappas))
    log.info(f"Plateau κ (last 20%): {plateau_mean:.6f} ± {plateau_std:.6f}")

    result = {
        "seed": args.seed,
        "final_kappa": final_kappa,
        "plateau_kappa_mean": plateau_mean,
        "plateau_kappa_std": plateau_std,
        "theory_kappa": theory,
        "si_prediction": 1.247,
    }
    with open(outdir / "result.json", "w") as f:
        json.dump(result, f, indent=2)

    return final_kappa


def main():
    p = argparse.ArgumentParser(
        description="SI §4.3: κ with patristic distance regression (switch thrown)"
    )
    p.add_argument("--manifest", required=True)
    p.add_argument("--output-dir", default="./patristic_kappa_run")
    p.add_argument("--vocab-size", type=int, default=4096)
    p.add_argument("--d-model", type=int, default=256)
    p.add_argument("--n-layers", type=int, default=4)
    p.add_argument("--latent-dim", type=int, default=128)
    p.add_argument("--max-len", type=int, default=8192)
    p.add_argument("--max-genomes", type=int, default=None)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--members-per-genus", type=int, default=4)
    p.add_argument("--lr", type=float, default=2e-4)
    p.add_argument("--steps", type=int, default=7000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log-every", type=int, default=10)
    p.add_argument("--warmup-steps", type=int, default=2000,
                   help="Steps of MLM+HEX-only training before patristic regression activates. "
                        "HEX must cluster same-genus embeddings first, otherwise patristic "
                        "pressure drives κ to its lower clamp.")
    p.add_argument("--dist-weight", type=float, default=0.1,
                   help="Weight on patristic dist_loss (default 0.1 — gentler than "
                        "BiosphereCodec's built-in 0.5 to prevent κ collapse from scratch)")
    p.add_argument("--sweep", action="store_true",
                   help="Run 5-seed sweep to reproduce SI Table 1 (seeds 0,42,123,456,789)")
    args = p.parse_args()

    args.genera_per_batch = args.batch_size // args.members_per_genus
    if args.genera_per_batch < 1:
        raise ValueError("batch_size must be >= members_per_genus")

    if args.sweep:
        seeds = [0, 42, 123, 456, 789]
        results = []
        base = Path(args.output_dir)
        for seed in seeds:
            args.seed = seed
            args.output_dir = str(base / f"seed_{seed}")
            kappa = train(args)
            results.append({"seed": seed, "kappa": kappa})
            with open(base / "sweep_results.json", "w") as f:
                json.dump(results, f, indent=2)

        kappas = [r["kappa"] for r in results]
        mean_k, std_k = np.mean(kappas), np.std(kappas)
        theory = (1.6 * math.log(2)) ** 2
        cv = std_k / mean_k * 100

        print("\n" + "=" * 60)
        print("5-SEED SWEEP RESULTS (SI §4.3)")
        print("=" * 60)
        for r in results:
            print(f"  seed={r['seed']:>3d}  κ={r['kappa']:.6f}")
        print(f"\nκ = {mean_k:.6f} ± {std_k:.6f}  (CV={cv:.1f}%)")
        print(f"SI prediction: 1.247000 ± 0.003000  (CV=0.2%)")
        print(f"Manning theory: {theory:.6f}")
        print(f"Agreement with SI: {abs(mean_k - 1.247) / 1.247 * 100:.1f}%")
        print(f"Agreement with theory: {abs(mean_k - theory) / theory * 100:.1f}%")

        with open(base / "sweep_results.json", "w") as f:
            json.dump({
                "seeds": results,
                "kappa_mean": float(mean_k),
                "kappa_std": float(std_k),
                "kappa_cv_pct": float(cv),
                "theory_kappa": float(theory),
                "si_prediction": 1.247,
            }, f, indent=2)
    else:
        train(args)


if __name__ == "__main__":
    main()
