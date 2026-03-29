#!/usr/bin/env python3
"""
Train BiosphereCodec from scratch → measure κ convergence
==========================================================

Phase 1: Pretrain BiosphereCodec on tokenized genomes (MLM + CLM + InfoNCE)
  - Encoder learns to compress genomic sequences (h enters via compression)
  - PoincareMapping projects to H^d with learnable curvature c
  - InfoNCE via tax_ids drives geometric organization

Phase 2: κ convergence test (5 seeds × multiple c_inits)
  - Freeze encoder, make c learnable from different starting points
  - Fine-tune with InfoNCE only
  - If κ converges to same value regardless of init → canonical constant

Usage:
  # Phase 1: Pretrain
  python train_codec_kappa.py --manifest /path/to/manifest.csv --phase pretrain

  # Phase 2: Measure κ (uses checkpoint from phase 1)
  python train_codec_kappa.py --manifest /path/to/manifest.csv --phase measure \
    --checkpoint pretrain_output/best.pt

  # Both phases
  python train_codec_kappa.py --manifest /path/to/manifest.csv --phase both
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Sampler

# ── Import BiosphereCodec ────────────────────────────────────────────────
# Add the model directory to path so we can import it
SCRIPT_DIR = Path(__file__).resolve().parent
CODEC_CANDIDATES = [
    SCRIPT_DIR.parent.parent / "model",  # ../../model relative to this script
    Path("/fast/sentrybio/scripts/Biosphere_codec/src/models"),
]

for p in CODEC_CANDIDATES:
    if (p / "BiosphereCodec.py").exists():
        sys.path.insert(0, str(p))
        break

from BiosphereCodec import BiosphereCodec


# ── Manual Poincaré distance (for κ measurement phase) ──────────────────

def manual_poincare_distance(u, v, c, eps=1e-7):
    diff_sq = torch.sum((u - v) ** 2, dim=-1)
    u_sq = torch.sum(u ** 2, dim=-1)
    v_sq = torch.sum(v ** 2, dim=-1)
    denom = (1.0 - c * u_sq).clamp(min=eps) * (1.0 - c * v_sq).clamp(min=eps)
    arg = 1.0 + 2.0 * c * diff_sq / denom
    return (1.0 / torch.sqrt(c + eps)) * torch.acosh(arg.clamp(min=1.0 + eps))


def patch_poincare_mapping(hyper_module):
    """Replace geoopt methods with manual Poincaré math for learnable c."""
    import types

    def patched_forward(self, x):
        z_euc = torch.tanh(self.lin(x))
        r_max = 0.9 / torch.sqrt(self.c.clamp(min=1e-7))
        z_euc = z_euc * r_max
        norms = torch.norm(z_euc, dim=-1, keepdim=True)
        scale = torch.clamp(r_max * 0.95 / (norms + 1e-7), max=1.0)
        return z_euc * scale

    def patched_dist_mat(self, z):
        B = z.shape[0]
        u = z.unsqueeze(1).expand(B, B, -1)
        v = z.unsqueeze(0).expand(B, B, -1)
        return manual_poincare_distance(u, v, self.c)

    hyper_module.forward = types.MethodType(patched_forward, hyper_module)
    hyper_module.dist_mat = types.MethodType(patched_dist_mat, hyper_module)
    hyper_module._manifold = None
    hyper_module._c_cached = None
    hyper_module._man = lambda: None


# ── Dataset ──────────────────────────────────────────────────────────────

class GenomeDataset(Dataset):
    """Load tokenized genomes with genus labels for InfoNCE."""

    def __init__(self, manifest_path, max_len=8192, max_genomes=None,
                 min_genus_count=2):
        import pandas as pd
        df = pd.read_csv(manifest_path, low_memory=False)
        df = df[df['tokenized_path'].notna() & df['genus'].notna()]
        df = df[df['tokenized_path'].apply(lambda p: os.path.exists(str(p)))]

        # Filter genera with enough members for InfoNCE
        gc = df['genus'].value_counts()
        valid = gc[gc >= min_genus_count].index
        df = df[df['genus'].isin(valid)]

        if max_genomes and len(df) > max_genomes:
            df = df.sample(max_genomes, random_state=42)

        self.paths = df['tokenized_path'].tolist()
        genera = df['genus'].tolist()
        self.genus_to_id = {g: i for i, g in enumerate(sorted(set(genera)))}
        self.labels = [self.genus_to_id[g] for g in genera]
        self.n_genera = len(self.genus_to_id)
        self.max_len = max_len

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        tokens = np.load(self.paths[idx]).astype(np.int64)
        if len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start:start + self.max_len]
        return torch.from_numpy(tokens), self.labels[idx]


def collate_genomes(batch):
    tokens_list, labels = zip(*batch)
    max_len = max(t.shape[0] for t in tokens_list)
    padded = torch.zeros(len(batch), max_len, dtype=torch.long)
    for i, t in enumerate(tokens_list):
        padded[i, :t.shape[0]] = t
    return padded, torch.tensor(labels, dtype=torch.long)


class GenusGroupedBatchSampler(Sampler):
    """Batch sampler that packs K genera × M members per batch.

    Guarantees every batch contains positive pairs for InfoNCE.
    With batch_size=32, members_per_genus=4: 8 genera × 4 = 32 samples,
    giving 8 × C(4,2) = 48 positive pairs per batch.
    """

    def __init__(self, labels: List[int], batch_size: int = 32,
                 members_per_genus: int = 4, seed: int = 42):
        self.batch_size = batch_size
        self.members_per_genus = members_per_genus
        self.genera_per_batch = batch_size // members_per_genus
        self.rng = np.random.RandomState(seed)

        # Build genus → sample indices mapping
        genus_indices: Dict[int, List[int]] = defaultdict(list)
        for idx, label in enumerate(labels):
            genus_indices[label].append(idx)

        # Only keep genera with enough members
        self.genus_indices = {
            g: idxs for g, idxs in genus_indices.items()
            if len(idxs) >= members_per_genus
        }
        self.genus_ids = list(self.genus_indices.keys())
        total_samples = sum(len(v) for v in self.genus_indices.values())
        self._len = max(1, total_samples // batch_size)

    def __iter__(self):
        genus_order = self.genus_ids.copy()
        self.rng.shuffle(genus_order)

        pools = {}
        for g in genus_order:
            idxs = self.genus_indices[g].copy()
            self.rng.shuffle(idxs)
            pools[g] = idxs

        g_ptr = 0
        while g_ptr + self.genera_per_batch <= len(genus_order):
            batch = []
            for _ in range(self.genera_per_batch):
                g = genus_order[g_ptr % len(genus_order)]
                g_ptr += 1
                pool = pools[g]
                for _ in range(self.members_per_genus):
                    if not pool:
                        pool = self.genus_indices[g].copy()
                        self.rng.shuffle(pool)
                        pools[g] = pool
                    batch.append(pool.pop())
            yield batch

    def __len__(self):
        return self._len


# ── Phase 1: Pretrain ────────────────────────────────────────────────────

def pretrain(args, log):
    """Train BiosphereCodec from scratch with MLM + CLM + InfoNCE."""
    device = torch.device(args.device)
    torch.manual_seed(42)
    random.seed(42)

    model = BiosphereCodec(
        vocab=args.vocab_size,
        d_model=args.d_model,
        n_layers=args.n_layers,
        max_len=args.max_len,
        latent_dim=args.latent_dim,
    )

    # Patch geoopt out BEFORE moving to device — avoids c corruption + NaN
    patch_poincare_mapping(model.hyper)
    model.hyper.c.data.fill_(1.0)  # reset c to 1.0 (geoopt may have corrupted it)
    model.hyper.c.requires_grad = True
    log(f"Patched PoincareMapping: manual Poincaré distance (c={model.hyper.c.item():.4f})")

    model = model.to(device)

    n_params = sum(p.numel() for p in model.parameters())
    log(f"BiosphereCodec: {n_params:,} params")
    log(f"  d_model={args.d_model}, n_layers={args.n_layers}, "
        f"latent_dim={args.latent_dim}")

    dataset = GenomeDataset(
        args.manifest, max_len=args.max_len,
        max_genomes=args.max_genomes,
        min_genus_count=args.min_genus_count,
    )
    log(f"Dataset: {len(dataset)} genomes, {dataset.n_genera} genera")

    sampler = GenusGroupedBatchSampler(
        dataset.labels, batch_size=args.batch_size,
        members_per_genus=args.members_per_genus, seed=42,
    )
    loader = DataLoader(
        dataset, batch_sampler=sampler,
        collate_fn=collate_genomes, num_workers=0,
    )
    n_genera_in_sampler = len(sampler.genus_ids)
    log(f"Grouped sampler: {sampler.genera_per_batch} genera × "
        f"{args.members_per_genus} members = {args.batch_size}/batch, "
        f"{n_genera_in_sampler} genera eligible")

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.pretrain_epochs,
    )

    best_loss = float('inf')
    outdir = Path(args.output_dir) / "pretrain"
    outdir.mkdir(parents=True, exist_ok=True)

    history = []
    t0 = time.time()

    for epoch in range(args.pretrain_epochs):
        model.train()
        epoch_loss = 0
        epoch_logs = defaultdict(float)
        n_batches = 0

        for tokens, labels in loader:
            tokens = tokens.to(device)
            labels = labels.to(device)

            loss, logs_dict = model(tokens, tax_ids=labels)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_loss += loss.item()
            for k, v in logs_dict.items():
                epoch_logs[k] += v
            n_batches += 1

        scheduler.step()

        avg_loss = epoch_loss / max(n_batches, 1)
        avg_logs = {k: v / max(n_batches, 1) for k, v in epoch_logs.items()}
        c_val = model.hyper.c.item()

        record = {
            'epoch': epoch, 'loss': avg_loss, 'c': c_val,
            **avg_logs, 'time': time.time() - t0,
        }
        history.append(record)

        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'loss': avg_loss,
                'c': c_val,
            }, outdir / "best.pt")

        if epoch % max(1, args.pretrain_epochs // 20) == 0 or epoch == args.pretrain_epochs - 1:
            log(f"  Ep {epoch:3d}/{args.pretrain_epochs}: loss={avg_loss:.4f} "
                f"mlm={avg_logs.get('mlm',0):.4f} dec={avg_logs.get('dec',0):.4f} "
                f"hex={avg_logs.get('hex',0):.4f} κ={c_val:.6f} "
                f"t={time.time()-t0:.0f}s")

    # Save final
    torch.save({
        'epoch': args.pretrain_epochs - 1,
        'model_state_dict': model.state_dict(),
        'loss': avg_loss,
        'c': model.hyper.c.item(),
    }, outdir / "final.pt")

    with open(outdir / "history.json", 'w') as f:
        json.dump(history, f, indent=2)

    log(f"\nPretrain complete. Best loss: {best_loss:.4f}")
    log(f"Checkpoint: {outdir / 'best.pt'}")
    log(f"Final κ (from pretraining): {model.hyper.c.item():.6f}")

    return str(outdir / "best.pt")


# ── Phase 2: κ Convergence Measurement ───────────────────────────────────

def measure_kappa(args, checkpoint_path, log):
    """Freeze encoder, sweep c_init × seeds, measure convergence."""
    device = torch.device(args.device)

    c_inits = [float(x) for x in args.c_inits.split(',')]
    seeds = [int(x) for x in args.seeds.split(',')]

    log(f"κ measurement: {len(seeds)} seeds × {len(c_inits)} inits "
        f"= {len(seeds)*len(c_inits)} runs")
    log(f"Checkpoint: {checkpoint_path}")
    log(f"c_inits: {c_inits}")
    log(f"seeds: {seeds}")
    log(f"Encoder: {'FROZEN' if args.freeze_encoder else 'fine-tuning'}")

    dataset = GenomeDataset(
        args.manifest, max_len=args.max_len,
        max_genomes=args.max_genomes,
        min_genus_count=args.min_genus_count,
    )
    log(f"Dataset: {len(dataset)} genomes, {dataset.n_genera} genera")

    outdir = Path(args.output_dir) / "kappa_measurement"
    outdir.mkdir(parents=True, exist_ok=True)

    all_results = []

    for seed in seeds:
        for c_init in c_inits:
            torch.manual_seed(seed)
            random.seed(seed)

            # Fresh model load each run
            model = BiosphereCodec(
                vocab=args.vocab_size,
                d_model=args.d_model,
                n_layers=args.n_layers,
                max_len=args.max_len,
                latent_dim=args.latent_dim,
            )

            ckpt = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
            state = ckpt.get('model_state_dict', ckpt)
            model.load_state_dict(state, strict=False)

            # Read pretrained c
            c_pretrained = model.hyper.c.item()

            # Patch to manual Poincaré (bypass geoopt)
            patch_poincare_mapping(model.hyper)

            # Set c to init value, make learnable
            model.hyper.c.data.fill_(c_init)
            model.hyper.c.requires_grad = True

            model = model.to(device)

            # Freeze encoder if requested
            if args.freeze_encoder:
                for name, p in model.named_parameters():
                    if name != 'hyper.c':
                        p.requires_grad = False

            # Optimizer: only c (and optionally encoder)
            param_groups = [
                {'params': [model.hyper.c], 'lr': args.c_lr, 'weight_decay': 0},
            ]
            if not args.freeze_encoder:
                encoder_params = [p for n, p in model.named_parameters()
                                  if n != 'hyper.c' and p.requires_grad]
                if encoder_params:
                    param_groups.append({
                        'params': encoder_params,
                        'lr': args.lr * 0.1,  # lower lr for fine-tuning
                        'weight_decay': 1e-4,
                    })

            optimizer = torch.optim.Adam(param_groups)

            sampler = GenusGroupedBatchSampler(
                dataset.labels, batch_size=args.batch_size,
                members_per_genus=args.members_per_genus, seed=seed,
            )
            loader = DataLoader(
                dataset, batch_sampler=sampler,
                collate_fn=collate_genomes, num_workers=0,
            )

            log(f"\n  seed={seed}, c_init={c_init:.2f} (pretrained={c_pretrained:.4f})")

            run_history = []
            t0 = time.time()

            for epoch in range(args.measure_epochs):
                model.train()
                epoch_loss = 0
                epoch_hex = 0
                n_batches = 0

                for tokens, labels in loader:
                    tokens = tokens.to(device)
                    labels = labels.to(device)

                    loss, logs_dict = model(tokens, tax_ids=labels)

                    optimizer.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    optimizer.step()

                    # Keep c positive
                    with torch.no_grad():
                        model.hyper.c.data.clamp_(min=0.01)

                    epoch_loss += loss.item()
                    epoch_hex += logs_dict.get('hex', 0)
                    n_batches += 1

                c_val = model.hyper.c.item()
                c_grad = model.hyper.c.grad.item() if model.hyper.c.grad is not None else 0
                avg_loss = epoch_loss / max(n_batches, 1)
                avg_hex = epoch_hex / max(n_batches, 1)

                run_history.append({
                    'epoch': epoch, 'c': c_val, 'c_grad': c_grad,
                    'loss': avg_loss, 'hex': avg_hex,
                })

                if epoch % max(1, args.measure_epochs // 10) == 0 or epoch == args.measure_epochs - 1:
                    log(f"    Ep {epoch:3d}: κ={c_val:.6f} ∇κ={c_grad:+.2e} "
                        f"loss={avg_loss:.4f} hex={avg_hex:.4f}")

            final_c = model.hyper.c.item()
            all_results.append({
                'seed': seed, 'c_init': c_init, 'c_final': final_c,
                'c_pretrained': c_pretrained,
                'loss_final': run_history[-1]['loss'],
                'hex_final': run_history[-1]['hex'],
            })
            log(f"    → Final κ = {final_c:.6f}")

            # Save run history
            with open(outdir / f"s{seed}_c{c_init:.2f}_history.json", 'w') as f:
                json.dump(run_history, f, indent=2)

            # Save incremental results
            with open(outdir / "results.json", 'w') as f:
                json.dump(all_results, f, indent=2)

    # ── Summary ──
    log("\n" + "=" * 60)
    log("κ CONVERGENCE SUMMARY")
    log("=" * 60)

    c_finals = np.array([r['c_final'] for r in all_results])
    mean_k = c_finals.mean()
    std_k = c_finals.std()
    cv = std_k / mean_k * 100 if mean_k > 0 else 999

    kappa_theory = (1.6 * math.log(2)) ** 2

    log(f"  {'seed':>6s}  {'c_init':>8s}  {'c_final':>10s}")
    log(f"  {'─'*6}  {'─'*8}  {'─'*10}")
    for r in all_results:
        log(f"  {r['seed']:>6d}  {r['c_init']:>8.2f}  {r['c_final']:>10.6f}")

    log(f"\n  κ = {mean_k:.6f} ± {std_k:.6f} (CV = {cv:.1f}%)")
    log(f"  Theory: κ = (1.6·ln2)² = {kappa_theory:.6f}")
    log(f"  Agreement: {abs(mean_k - kappa_theory) / kappa_theory * 100:.1f}%")

    # Convergence by init
    log(f"\n  By init_κ:")
    for c_init in c_inits:
        group = [r['c_final'] for r in all_results if r['c_init'] == c_init]
        if group:
            log(f"    init={c_init:.2f}: {np.mean(group):.6f} ± {np.std(group):.6f}")

    with open(outdir / "summary.json", 'w') as f:
        json.dump({
            'method': 'codec_kappa_finetune',
            'latent_dim': args.latent_dim,
            'freeze_encoder': args.freeze_encoder,
            'pretrain_epochs': args.pretrain_epochs,
            'measure_epochs': args.measure_epochs,
            'n_genomes': len(dataset),
            'n_genera': dataset.n_genera,
            'kappa_mean': float(mean_k),
            'kappa_std': float(std_k),
            'kappa_cv': float(cv),
            'kappa_theory': float(kappa_theory),
            'agreement_pct': float(abs(mean_k - kappa_theory) / kappa_theory * 100),
            'results': all_results,
        }, f, indent=2)


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Train BiosphereCodec → Measure κ convergence")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="codec_kappa_output")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--phase", choices=["pretrain", "measure", "both"], default="both")
    parser.add_argument("--checkpoint", type=str, default=None,
                        help="Checkpoint for measure phase (skip pretrain)")

    # Architecture
    parser.add_argument("--vocab-size", type=int, default=4096)
    parser.add_argument("--d-model", type=int, default=256,
                        help="Model dimension (256 for speed, 512 for quality)")
    parser.add_argument("--n-layers", type=int, default=4,
                        help="Encoder layers (4 for speed, 6 for quality)")
    parser.add_argument("--latent-dim", type=int, default=2,
                        help="Poincaré ball dimension (2 for H²)")
    parser.add_argument("--max-len", type=int, default=4096,
                        help="Max sequence length")

    # Data
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--min-genus-count", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=32,
                        help="Batch size (genera_per_batch × members_per_genus)")
    parser.add_argument("--members-per-genus", type=int, default=4,
                        help="Members per genus in grouped sampler")

    # Pretrain
    parser.add_argument("--pretrain-epochs", type=int, default=30)
    parser.add_argument("--lr", type=float, default=3e-4)

    # κ measurement
    parser.add_argument("--measure-epochs", type=int, default=20)
    parser.add_argument("--c-lr", type=float, default=1e-2)
    parser.add_argument("--c-inits", type=str, default="0.5,1.0,1.25,2.0,5.0")
    parser.add_argument("--seeds", type=str, default="42,59,76,93,110")
    parser.add_argument("--freeze-encoder", action="store_true", default=True)
    parser.add_argument("--no-freeze-encoder", dest="freeze_encoder", action="store_false")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "train_codec_kappa.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("BIOSPHERE CODEC — TRAIN + κ CONVERGENCE TEST")
    log("=" * 70)
    log(f"Phase: {args.phase}")
    log(f"Architecture: d_model={args.d_model}, n_layers={args.n_layers}, "
        f"latent_dim={args.latent_dim}")
    log(f"Device: {args.device}")
    log()

    checkpoint_path = args.checkpoint

    try:
        if args.phase in ("pretrain", "both"):
            log("=" * 60)
            log("PHASE 1: PRETRAIN")
            log("=" * 60)
            checkpoint_path = pretrain(args, log)
            log()

        if args.phase in ("measure", "both"):
            if checkpoint_path is None:
                log("ERROR: --checkpoint required for measure phase")
                return
            log("=" * 60)
            log("PHASE 2: κ CONVERGENCE MEASUREMENT")
            log("=" * 60)
            measure_kappa(args, checkpoint_path, log)
    except Exception:
        import traceback
        tb = traceback.format_exc()
        log(f"\nFATAL ERROR:\n{tb}")
        raise
    finally:
        log_file.close()


if __name__ == "__main__":
    main()
