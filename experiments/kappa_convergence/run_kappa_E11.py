#!/usr/bin/env python3
"""
κ Convergence via E11 Architecture — Learnable Curvature on 38K Genomes
========================================================================

This is the E11 coordinate experiment with ONE change: κ is learnable.

Architecture: E11's PoincareEncoder2D (40K params) with exp_map using
learnable c instead of frozen KAPPA=1.25.

Loss: Same 5-term loss as E11 (quartet + domain + genus + repel + radial),
all of which compute Poincare distances that depend on c.

Data: Pre-tokenized .npy files via manifest_local.csv (38K genomes, 268 anchors).

Hypothesis: If κ=1.25 is a property of the tree of life, it should emerge
from the data regardless of initialization.

Usage:
  # Single test
  python run_kappa_E11.py --manifest /path/to/manifest.csv --init-kappa 0.5

  # Full sweep: 5 seeds × 3 initializations
  python run_kappa_E11.py --manifest /path/to/manifest.csv --sweep
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

# ── Constants (same as E11 EXCEPT kappa is now learnable) ─────────────────

VOCAB_SIZE = 4096
FEAT_DIM = 8
MAX_SEQ_LEN = 512     # tokenized .npy files are 512 tokens
LATENT_DIM = 2

TAXONOMY_RANKS = ["genus", "family", "order", "class", "phylum", "domain"]

DOMAIN_ANGLES = {
    "Bacteria":   0.0,
    "Archaea":    120.0,
    "Eukaryota":  240.0,
}

GENUS_ANCHORS = {
    "Escherichia":         0.0,
    "Methanocaldococcus":  120.0,
    "Saccharomyces":       240.0,
}


# ── Poincaré Ball Operations (LEARNABLE curvature) ───────────────────────

def poincare_distance(x: torch.Tensor, y: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Geodesic distance in the Poincaré ball with learnable curvature c."""
    c_val = torch.clamp(c, min=1e-4)
    sqrt_c = torch.sqrt(c_val)
    x_sq = (x * x).sum(-1, keepdim=True).clamp(max=1.0 / c_val.item() - 1e-5)
    y_sq = (y * y).sum(-1, keepdim=True).clamp(max=1.0 / c_val.item() - 1e-5)
    xy_diff_sq = ((x - y) ** 2).sum(-1, keepdim=True)
    num = 2.0 * c_val * xy_diff_sq
    denom = (1.0 - c_val * x_sq) * (1.0 - c_val * y_sq)
    arg = 1.0 + num / denom.clamp(min=1e-10)
    return ((1.0 / sqrt_c) * torch.acosh(arg.clamp(min=1.0 + 1e-7))).squeeze(-1)


def exp_map_zero(v: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Exponential map at the origin with learnable curvature."""
    c_val = torch.clamp(c, min=1e-4)
    sqrt_c = torch.sqrt(c_val)
    v_norm = v.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    return torch.tanh(sqrt_c * v_norm / 2.0) * v / (sqrt_c * v_norm)


def project_to_ball(x: torch.Tensor, c: torch.Tensor, eps: float = 1e-5) -> torch.Tensor:
    """Project onto open ball {x : c||x||^2 < 1}."""
    c_val = torch.clamp(c, min=1e-4)
    max_norm = (1.0 / torch.sqrt(c_val)) - eps
    norms = x.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    clamped = (max_norm / norms).clamp(max=1.0)
    return x * clamped


# ── Model ─────────────────────────────────────────────────────────────────

class PoincareEncoder2DLearnable(nn.Module):
    """E11 encoder with LEARNABLE curvature c as nn.Parameter."""

    CONV_DIM = 32

    def __init__(self, init_c: float = 1.25, vocab_size: int = VOCAB_SIZE,
                 feat_dim: int = FEAT_DIM, hidden: int = 64):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, feat_dim)
        self.conv = nn.Conv1d(feat_dim, self.CONV_DIM, kernel_size=5, padding=2)
        self.encoder = nn.Sequential(
            nn.Linear(self.CONV_DIM, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, LATENT_DIM),
        )
        # THE KEY: learnable curvature
        self.c = nn.Parameter(torch.tensor(init_c))

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        x = self.embed(tokens)
        x = x.transpose(1, 2)
        x = F.gelu(self.conv(x))
        x = x.mean(dim=2)
        v = self.encoder(x)
        z = exp_map_zero(v, self.c)
        return project_to_ball(z, self.c)


# ── Dataset (reused from E11) ────────────────────────────────────────────

class AnchorDataset(Dataset):
    def __init__(self, anchors: List[Dict], max_len: int = MAX_SEQ_LEN,
                 augment: bool = True):
        self.anchors = anchors
        self.max_len = max_len
        self.augment = augment

    def __len__(self) -> int:
        return len(self.anchors)

    def __getitem__(self, idx: int) -> Dict:
        a = self.anchors[idx]
        path = a["tokenized_path"]
        tokens = torch.from_numpy(np.load(path)).long()

        if self.augment and len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start:start + self.max_len]
        elif len(tokens) > self.max_len:
            tokens = tokens[:self.max_len]

        return {
            "tokens": tokens,
            "idx": idx,
            "domain": a.get("domain", "Unknown"),
            "phylum": a.get("phylum", "Unknown"),
            "class": a.get("class", "Unknown"),
            "order": a.get("order", "Unknown"),
            "family": a.get("family", "Unknown"),
            "genus": a.get("genus", "Unknown"),
            "genome_size": float(a.get("genome_size", 0) or 0),
        }


def collate_anchors(batch: List[Dict]) -> Dict:
    max_len = max(b["tokens"].shape[0] for b in batch)
    padded = torch.zeros(len(batch), max_len, dtype=torch.long)
    for i, b in enumerate(batch):
        L = b["tokens"].shape[0]
        padded[i, :L] = b["tokens"]
    return {
        "tokens": padded,
        "idx": torch.tensor([b["idx"] for b in batch]),
        "domains": [b["domain"] for b in batch],
        "phyla": [b["phylum"] for b in batch],
        "classes": [b["class"] for b in batch],
        "orders": [b["order"] for b in batch],
        "families": [b["family"] for b in batch],
        "genera": [b["genus"] for b in batch],
        "genome_sizes": torch.tensor([b["genome_size"] for b in batch], dtype=torch.float32),
    }


# ── Anchor Selection (from E11) ──────────────────────────────────────────

def select_anchors(manifest_path: str, n_target: int = 250,
                   seed: int = 42) -> List[Dict]:
    rng = random.Random(seed)
    all_rows = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            tp = row.get("tokenized_path", "")
            if not tp or not os.path.exists(tp):
                continue
            if row.get("domain") in ("Viruses", ""):
                continue
            all_rows.append(row)

    tree = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in all_rows:
        d = row.get("domain", "Unknown")
        p = row.get("phylum", "Unknown")
        f = row.get("family", "Unknown")
        tree[d][p][f].append(row)

    domains = sorted(tree.keys())
    domain_budget = {}
    for d in domains:
        n_available = sum(len(fam) for phy in tree[d].values() for fam in phy.values())
        domain_budget[d] = max(30, int(n_target * n_available / len(all_rows)))

    anchors = []
    for domain in domains:
        budget = domain_budget[domain]
        phyla = sorted(tree[domain].keys())
        n_per_phylum = max(2, budget // max(len(phyla), 1))
        for phylum in phyla:
            families = sorted(tree[domain][phylum].keys())
            n_per_family = max(1, n_per_phylum // max(len(families), 1))
            for family in families:
                candidates = tree[domain][phylum][family]
                candidates.sort(key=lambda r: float(r.get("genome_size", 0) or 0))
                if len(candidates) <= n_per_family:
                    anchors.extend(candidates)
                else:
                    indices = np.linspace(0, len(candidates) - 1, n_per_family, dtype=int)
                    anchors.extend([candidates[i] for i in indices])

    # Ensure genus anchors included
    anchor_accessions = {a.get("accession") for a in anchors}
    for genus in GENUS_ANCHORS:
        found = any(a.get("genus") == genus for a in anchors)
        if not found:
            for row in all_rows:
                if row.get("genus") == genus and row.get("accession") not in anchor_accessions:
                    anchors.append(row)
                    break

    if len(anchors) > n_target * 2:
        rng.shuffle(anchors)
        trimmed = []
        by_domain = defaultdict(list)
        for a in anchors:
            by_domain[a.get("domain", "Unknown")].append(a)
        for d in domains:
            pool = by_domain[d]
            take = max(30, int(n_target * len(pool) / len(anchors)))
            trimmed.extend(pool[:take])
        anchors = trimmed

    return anchors


# ── Taxonomy Distance ─────────────────────────────────────────────────────

def batch_taxonomy_distance(batch: Dict, i: int, j: int) -> int:
    for level, rank_key in enumerate(
        ["genera", "families", "orders", "classes", "phyla", "domains"]
    ):
        va = batch[rank_key][i]
        vb = batch[rank_key][j]
        if va == vb and va not in ("Unknown", ""):
            return level
    return 6


# ── Loss Functions (same as E11 but using learnable c) ────────────────────

def quartet_loss(coords: torch.Tensor, batch: Dict, c: torch.Tensor,
                 n_quartets: int = 200, margin: float = 0.5) -> torch.Tensor:
    B = coords.shape[0]
    if B < 4:
        return torch.tensor(0.0, device=coords.device)

    indices = list(range(B))
    n_quartets = min(n_quartets, B * 3)
    loss = torch.tensor(0.0, device=coords.device)
    count = 0

    for _ in range(n_quartets):
        q = random.sample(indices, 4)
        splits = [
            (q[0], q[1], q[2], q[3]),
            (q[0], q[2], q[1], q[3]),
            (q[0], q[3], q[1], q[2]),
        ]
        best = min(splits, key=lambda s: (
            batch_taxonomy_distance(batch, s[0], s[1])
            + batch_taxonomy_distance(batch, s[2], s[3])
        ))
        i, j, k, l = best
        d_ij = poincare_distance(coords[i:i+1], coords[j:j+1], c)
        d_kl = poincare_distance(coords[k:k+1], coords[l:l+1], c)
        d_ik = poincare_distance(coords[i:i+1], coords[k:k+1], c)
        d_jl = poincare_distance(coords[j:j+1], coords[l:l+1], c)
        violation = (d_ij + d_kl) - (d_ik + d_jl) + margin
        loss = loss + F.relu(violation).squeeze()
        count += 1

    return loss / max(count, 1)


def domain_angular_loss(coords: torch.Tensor, domains: List[str]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    target_dirs = {}
    for domain, angle_deg in DOMAIN_ANGLES.items():
        rad = math.radians(angle_deg)
        target_dirs[domain] = torch.tensor([math.cos(rad), math.sin(rad)], device=coords.device)
    for i, domain in enumerate(domains):
        if domain in target_dirs:
            coord_norm = coords[i] / coords[i].norm().clamp(min=1e-8)
            cos_sim = (coord_norm * target_dirs[domain]).sum()
            loss = loss + (1.0 - cos_sim)
            count += 1
    return loss / max(count, 1)


def genus_anchor_loss(coords: torch.Tensor, genera: List[str]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for i, genus in enumerate(genera):
        if genus in GENUS_ANCHORS:
            target_rad = math.radians(GENUS_ANCHORS[genus])
            target_dir = torch.tensor([math.cos(target_rad), math.sin(target_rad)], device=coords.device)
            coord_norm = coords[i] / coords[i].norm().clamp(min=1e-8)
            cos_sim = (coord_norm * target_dir).sum()
            loss = loss + (1.0 - cos_sim)
            count += 1
    return loss / max(count, 1)


def angular_repulsion_loss(coords: torch.Tensor, families: List[str],
                           domains: List[str], n_pairs: int = 300) -> torch.Tensor:
    B = coords.shape[0]
    if B < 2:
        return torch.tensor(0.0, device=coords.device)
    norms = coords.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    dirs = coords / norms
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    indices = list(range(B))
    for _ in range(n_pairs):
        i, j = random.sample(indices, 2)
        if families[i] == families[j] and families[i] not in ("Unknown", ""):
            continue
        cos_sim = (dirs[i] * dirs[j]).sum()
        if domains[i] != domains[j]:
            margin = 0.5
        else:
            margin = 0.966
        loss = loss + F.relu(cos_sim - margin)
        count += 1
    return loss / max(count, 1)


def radial_ordering_loss(coords: torch.Tensor, genome_sizes: torch.Tensor,
                         n_pairs: int = 300, margin: float = 0.02) -> torch.Tensor:
    r = coords.norm(dim=-1)
    B = r.shape[0]
    if B < 2:
        return torch.tensor(0.0, device=coords.device)
    valid = (genome_sizes > 0).nonzero(as_tuple=True)[0].tolist()
    if len(valid) < 2:
        return torch.tensor(0.0, device=coords.device)
    n_pairs = min(n_pairs, len(valid) * 3)
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for _ in range(n_pairs):
        i, j = random.sample(valid, 2)
        si, sj = genome_sizes[i].item(), genome_sizes[j].item()
        if si > sj * 1.5:
            loss = loss + F.relu(r[j] - r[i] + margin)
            count += 1
        elif sj > si * 1.5:
            loss = loss + F.relu(r[i] - r[j] + margin)
            count += 1
    return loss / max(count, 1)


# ── Training ──────────────────────────────────────────────────────────────

def train_one_run(
    anchors: List[Dict],
    seed: int,
    init_c: float,
    n_epochs: int,
    lr: float,
    batch_size: int,
    device: str,
    output_dir: str,
    log,
) -> Tuple[float, List[Dict]]:
    """Train one model with learnable κ. Return final κ and history."""
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    model = PoincareEncoder2DLearnable(init_c=init_c).to(device)

    # Separate optimizer groups: κ gets no weight decay
    other_params = [p for n, p in model.named_parameters() if n != 'c']
    optimizer = torch.optim.AdamW([
        {'params': other_params, 'lr': lr, 'weight_decay': 1e-4},
        {'params': [model.c], 'lr': lr, 'weight_decay': 0.0},
    ])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    dataset = AnchorDataset(anchors, augment=True)
    loader = DataLoader(
        dataset, batch_size=min(len(anchors), batch_size),
        shuffle=True, collate_fn=collate_anchors, num_workers=0,
    )

    n_params = sum(p.numel() for p in model.parameters())
    log(f"  Seed {seed}, init_κ={init_c:.2f}: {n_params:,} params, c.requires_grad={model.c.requires_grad}")

    kappa_history = []

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for batch in loader:
            tokens = batch["tokens"].to(device)
            coords = model(tokens)

            lq = quartet_loss(coords, batch, model.c)
            la = domain_angular_loss(coords, batch["domains"])
            lg = genus_anchor_loss(coords, batch["genera"])
            lrep = angular_repulsion_loss(coords, batch["families"], batch["domains"])
            lr_loss = radial_ordering_loss(coords, batch["genome_sizes"].to(device))

            loss = lq + 2.0 * la + 5.0 * lg + 0.5 * lrep + 0.3 * lr_loss

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()

        kappa_val = model.c.item()
        kappa_grad = model.c.grad.item() if model.c.grad is not None else 0.0
        avg_loss = epoch_loss / max(n_batches, 1)

        kappa_history.append({
            'epoch': epoch + 1,
            'kappa': kappa_val,
            'kappa_grad': kappa_grad,
            'loss': avg_loss,
        })

        if (epoch + 1) % 10 == 0 or epoch == 0 or epoch == n_epochs - 1:
            log(f"    Epoch {epoch+1:4d}: κ={kappa_val:.6f}  ∇κ={kappa_grad:+.2e}  loss={avg_loss:.4f}")

    final_kappa = model.c.item()

    # Save checkpoint and history
    run_dir = Path(output_dir) / f"seed_{seed}_init_{init_c:.2f}"
    run_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), run_dir / "model.pt")
    with open(run_dir / "kappa_history.json", 'w') as f:
        json.dump(kappa_history, f, indent=2)

    return final_kappa, kappa_history


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="κ Convergence via E11 Architecture")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="kappa_E11_results")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--n-anchors", type=int, default=250)
    parser.add_argument("--n-epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--init-kappa", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sweep", action="store_true",
                        help="Run 5 seeds × 3 initializations")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "kappa_E11.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("κ CONVERGENCE VIA E11 ARCHITECTURE (LEARNABLE CURVATURE)")
    log("=" * 70)
    log(f"Manifest: {args.manifest}")
    log(f"Device: {args.device}")
    log()

    # Select anchors
    log("--- Selecting anchor organisms ---")
    anchors = select_anchors(args.manifest, n_target=args.n_anchors)
    domain_counts = defaultdict(int)
    for a in anchors:
        domain_counts[a.get("domain", "Unknown")] += 1
    log(f"  Selected {len(anchors)} anchors: {dict(domain_counts)}")
    log()

    if args.sweep:
        seeds = [42, 59, 76, 93, 110]  # Same as E11
        init_kappas = [0.5, 1.0, 2.0]
        results = []

        for seed in seeds:
            for init_k in init_kappas:
                log(f"{'='*60}")
                log(f"  SWEEP: seed={seed}, init_κ={init_k}")
                log(f"{'='*60}")

                final_k, history = train_one_run(
                    anchors, seed, init_k, args.n_epochs, args.lr,
                    args.batch_size, args.device, args.output_dir, log,
                )
                results.append({
                    'seed': seed, 'init_kappa': init_k,
                    'final_kappa': final_k,
                })
                log(f"  → Final κ = {final_k:.6f}")
                log()

                # Save intermediate
                with open(os.path.join(args.output_dir, "sweep_results.json"), 'w') as f:
                    json.dump(results, f, indent=2)

        # Summary
        log("=" * 70)
        log("SWEEP COMPLETE")
        log("=" * 70)
        for r in results:
            log(f"  seed={r['seed']:>4d}  init_κ={r['init_kappa']:.1f}  →  final_κ={r['final_kappa']:.6f}")

        # Group by init_kappa
        for init_k in init_kappas:
            kappas = [r['final_kappa'] for r in results if r['init_kappa'] == init_k]
            log(f"\n  Init {init_k:.1f}: mean={np.mean(kappas):.6f} ± {np.std(kappas):.6f}")

        all_kappas = [r['final_kappa'] for r in results]
        log(f"\n  Overall: mean={np.mean(all_kappas):.6f} ± {np.std(all_kappas):.6f}")
        log(f"  Range: [{min(all_kappas):.6f}, {max(all_kappas):.6f}]")

    else:
        final_k, history = train_one_run(
            anchors, args.seed, args.init_kappa, args.n_epochs, args.lr,
            args.batch_size, args.device, args.output_dir, log,
        )
        log(f"\nFinal κ = {final_k:.6f} (init was {args.init_kappa})")

    log_file.close()


if __name__ == "__main__":
    main()
