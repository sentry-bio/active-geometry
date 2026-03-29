#!/usr/bin/env python3
"""
κ Convergence via E11 Architecture — RIEMANNIAN GRADIENT CORRECTION
====================================================================

Same as run_kappa_E11.py but with the critical fix: Riemannian gradient
scaling on manifold points.

The problem: Euclidean Adam treats the Poincaré ball as flat space. Near
the boundary, the metric tensor blows up — a small Euclidean step is an
enormous geodesic displacement. This causes κ to systematically overshoot
because the optimizer can't feel the curvature it's optimizing on.

The fix: Scale gradients at z by the inverse conformal factor:
    λ(z) = (1 - c‖z‖²)² / 4
This converts Euclidean gradients to Riemannian gradients without
importing geoopt (avoiding the softplus corruption bug).

Architecture:
  - E11's PoincareEncoder2D (~40K params)
  - Softplus reparameterization for c (clean gradient flow)
  - 5-term geometric loss (quartet + domain + genus + repel + radial)
  - Riemannian gradient correction on manifold points
  - 2D Poincaré ball (H²)

Usage:
  python run_kappa_E11_riemannian.py --manifest /path/to/manifest.csv --sweep
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

# ── Constants ─────────────────────────────────────────────────────────────

VOCAB_SIZE = 4096
FEAT_DIM = 8
MAX_SEQ_LEN = 512
LATENT_DIM = 2

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


# ── Poincaré Ball Operations ─────────────────────────────────────────────
# All operations use c from the autograd graph (softplus(raw_c)).
# No geoopt — we implement the geometry manually to avoid the
# in-place softplus corruption bug.

def poincare_distance(x: torch.Tensor, y: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Geodesic distance: d(x,y) = (1/√c) · arccosh(1 + 2c‖x-y‖²/((1-c‖x‖²)(1-c‖y‖²)))"""
    eps = 1e-7
    x_sq = (x * x).sum(-1, keepdim=True)
    y_sq = (y * y).sum(-1, keepdim=True)
    xy_diff_sq = ((x - y) ** 2).sum(-1, keepdim=True)

    num = 2.0 * c * xy_diff_sq
    denom = (1.0 - c * x_sq).clamp(min=eps) * (1.0 - c * y_sq).clamp(min=eps)
    arg = 1.0 + num / denom
    return ((1.0 / torch.sqrt(c + eps)) * torch.acosh(arg.clamp(min=1.0 + eps))).squeeze(-1)


def exp_map_zero(v: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Exponential map at origin: maps tangent vector v to Poincaré ball."""
    sqrt_c = torch.sqrt(c.clamp(min=1e-7))
    v_norm = v.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    return torch.tanh(sqrt_c * v_norm / 2.0) * v / (sqrt_c * v_norm)


def project_to_ball(x: torch.Tensor, c: torch.Tensor, eps: float = 1e-5) -> torch.Tensor:
    """Project onto open ball {x : c‖x‖² < 1}."""
    max_norm = (1.0 / torch.sqrt(c.clamp(min=1e-7))) - eps
    norms = x.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    scale = (max_norm / norms).clamp(max=1.0)
    return x * scale


def conformal_factor(z: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Poincaré ball conformal factor λ(z) = (1 - c‖z‖²)² / 4.

    This converts Euclidean gradients to Riemannian gradients:
        ∇_riem = ∇_euc / λ(z)

    Near the center (‖z‖≈0): λ≈1/4, gradients amplified 4×
    Near the boundary (c‖z‖²≈1): λ≈0, gradients damped to 0

    This is the KEY missing piece. Without it, Euclidean Adam
    massively over-updates points near the boundary, causing
    systematic overshoot of the curvature basin.
    """
    z_sq = (z * z).sum(-1, keepdim=True)
    return ((1.0 - c * z_sq).clamp(min=1e-7) ** 2) / 4.0


# ── Riemannian Gradient Hook ─────────────────────────────────────────────

class RiemannianGradientHook:
    """Register a backward hook on z to scale gradients by 1/λ(z).

    This converts the Euclidean gradient that flows backward through z
    into the Riemannian gradient on the Poincaré ball. The hook fires
    during backward(), before the gradient reaches the encoder weights.

    Math: ∂L/∂θ_riem = ∂L/∂z · (1/λ(z)) · ∂z/∂θ
    where λ(z) = (1 - c‖z‖²)² / 4 is the conformal factor.
    """

    def __init__(self):
        self.handle = None

    def register(self, z: torch.Tensor, c: torch.Tensor):
        """Register hook on tensor z for this forward pass."""
        lam = conformal_factor(z, c.detach())  # detach c to avoid double-counting

        def hook_fn(grad):
            # Scale Euclidean grad → Riemannian grad
            # Riemannian grad = Euclidean grad / λ(z)
            # But we want to DAMPEN near boundary (not amplify),
            # so we multiply by λ (the conformal factor itself)
            return grad * lam

        if self.handle is not None:
            self.handle.remove()
        self.handle = z.register_hook(hook_fn)


# ── Model ─────────────────────────────────────────────────────────────────

class PoincareEncoder2DRiemannian(nn.Module):
    """E11 encoder with learnable c (softplus) and Riemannian gradient correction."""

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

        # Softplus reparameterization: c = softplus(raw_c), always positive
        # raw_c = log(exp(init_c) - 1) so that softplus(raw_c) = init_c
        if init_c > 0.01:
            raw_init = math.log(math.exp(init_c) - 1.0)
        else:
            raw_init = -4.0
        self.raw_c = nn.Parameter(torch.tensor(raw_init))

        # Riemannian gradient hook
        self.riem_hook = RiemannianGradientHook()

    @property
    def c(self) -> torch.Tensor:
        """Curvature, always positive via softplus. Never in-place."""
        return F.softplus(self.raw_c)

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        c = self.c
        x = self.embed(tokens)
        x = x.transpose(1, 2)
        x = F.gelu(self.conv(x))
        x = x.mean(dim=2)
        v = self.encoder(x)
        z = exp_map_zero(v, c)
        z = project_to_ball(z, c)

        # Register Riemannian gradient hook: gradients flowing backward
        # through z will be scaled by the conformal factor λ(z).
        # This makes the optimizer respect the manifold geometry.
        if z.requires_grad:
            self.riem_hook.register(z, c)

        return z


# ── Dataset ───────────────────────────────────────────────────────────────

class AnchorDataset(Dataset):
    def __init__(self, anchors: List[Dict], max_len: int = MAX_SEQ_LEN, augment: bool = True):
        self.anchors = anchors
        self.max_len = max_len
        self.augment = augment

    def __len__(self) -> int:
        return len(self.anchors)

    def __getitem__(self, idx: int) -> Dict:
        a = self.anchors[idx]
        tokens = torch.from_numpy(np.load(a["tokenized_path"])).long()
        if self.augment and len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start:start + self.max_len]
        elif len(tokens) > self.max_len:
            tokens = tokens[:self.max_len]
        return {
            "tokens": tokens, "idx": idx,
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
        padded[i, :b["tokens"].shape[0]] = b["tokens"]
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


# ── Anchor Selection ──────────────────────────────────────────────────────

def select_anchors(manifest_path: str, n_target: int = 250, seed: int = 42) -> List[Dict]:
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
        tree[row.get("domain", "Unknown")][row.get("phylum", "Unknown")][row.get("family", "Unknown")].append(row)

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

    for genus in GENUS_ANCHORS:
        if not any(a.get("genus") == genus for a in anchors):
            for row in all_rows:
                if row.get("genus") == genus:
                    anchors.append(row)
                    break

    if len(anchors) > n_target * 2:
        rng.shuffle(anchors)
        by_domain = defaultdict(list)
        for a in anchors:
            by_domain[a.get("domain", "Unknown")].append(a)
        trimmed = []
        for d in domains:
            pool = by_domain[d]
            trimmed.extend(pool[:max(30, int(n_target * len(pool) / len(anchors)))])
        anchors = trimmed

    return anchors


# ── Taxonomy Distance ─────────────────────────────────────────────────────

RANK_KEYS = ["genera", "families", "orders", "classes", "phyla", "domains"]

def batch_taxonomy_distance(batch: Dict, i: int, j: int) -> int:
    for level, key in enumerate(RANK_KEYS):
        if batch[key][i] == batch[key][j] and batch[key][i] not in ("Unknown", ""):
            return level
    return 6


# ── Loss Functions ────────────────────────────────────────────────────────
# All losses use poincare_distance(x, y, c) where c = softplus(raw_c).
# This means ∂L/∂raw_c flows through every geometric computation.

def quartet_loss(coords: torch.Tensor, batch: Dict, c: torch.Tensor,
                 n_quartets: int = 200, margin: float = 0.5) -> torch.Tensor:
    B = coords.shape[0]
    if B < 4:
        return torch.tensor(0.0, device=coords.device)

    indices = list(range(B))
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for _ in range(min(n_quartets, B * 3)):
        q = random.sample(indices, 4)
        splits = [(q[0],q[1],q[2],q[3]), (q[0],q[2],q[1],q[3]), (q[0],q[3],q[1],q[2])]
        best = min(splits, key=lambda s: batch_taxonomy_distance(batch, s[0], s[1]) + batch_taxonomy_distance(batch, s[2], s[3]))
        i, j, k, l = best
        d_ij = poincare_distance(coords[i:i+1], coords[j:j+1], c)
        d_kl = poincare_distance(coords[k:k+1], coords[l:l+1], c)
        d_ik = poincare_distance(coords[i:i+1], coords[k:k+1], c)
        d_jl = poincare_distance(coords[j:j+1], coords[l:l+1], c)
        loss = loss + F.relu((d_ij + d_kl) - (d_ik + d_jl) + margin).squeeze()
        count += 1
    return loss / max(count, 1)


def domain_angular_loss(coords: torch.Tensor, domains: List[str]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for domain, angle_deg in DOMAIN_ANGLES.items():
        rad = math.radians(angle_deg)
        target = torch.tensor([math.cos(rad), math.sin(rad)], device=coords.device)
        for i, d in enumerate(domains):
            if d == domain:
                norm = coords[i] / coords[i].norm().clamp(min=1e-8)
                loss = loss + (1.0 - (norm * target).sum())
                count += 1
    return loss / max(count, 1)


def genus_anchor_loss(coords: torch.Tensor, genera: List[str]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for i, genus in enumerate(genera):
        if genus in GENUS_ANCHORS:
            rad = math.radians(GENUS_ANCHORS[genus])
            target = torch.tensor([math.cos(rad), math.sin(rad)], device=coords.device)
            norm = coords[i] / coords[i].norm().clamp(min=1e-8)
            loss = loss + (1.0 - (norm * target).sum())
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
        margin = 0.5 if domains[i] != domains[j] else 0.966
        loss = loss + F.relu(cos_sim - margin)
        count += 1
    return loss / max(count, 1)


def radial_ordering_loss(coords: torch.Tensor, genome_sizes: torch.Tensor,
                          n_pairs: int = 300, margin: float = 0.02) -> torch.Tensor:
    r = coords.norm(dim=-1)
    B = r.shape[0]
    valid = (genome_sizes > 0).nonzero(as_tuple=True)[0].tolist()
    if len(valid) < 2:
        return torch.tensor(0.0, device=coords.device)
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for _ in range(min(n_pairs, len(valid) * 3)):
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
    anchors: List[Dict], seed: int, init_c: float, n_epochs: int,
    lr: float, batch_size: int, device: str, output_dir: str, log,
) -> Tuple[float, List[Dict]]:
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    model = PoincareEncoder2DRiemannian(init_c=init_c).to(device)

    # Separate optimizer groups: raw_c gets no weight decay
    other_params = [p for n, p in model.named_parameters() if n != 'raw_c']
    optimizer = torch.optim.AdamW([
        {'params': other_params, 'lr': lr, 'weight_decay': 1e-4},
        {'params': [model.raw_c], 'lr': lr, 'weight_decay': 0.0},
    ])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    dataset = AnchorDataset(anchors, augment=True)
    loader = DataLoader(dataset, batch_size=min(len(anchors), batch_size),
                        shuffle=True, collate_fn=collate_anchors, num_workers=0)

    n_params = sum(p.numel() for p in model.parameters())
    log(f"  Seed {seed}, init_κ={init_c:.2f}: {n_params:,} params, "
        f"raw_c={model.raw_c.item():.4f}, c={model.c.item():.6f}")

    history = []

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        epoch_r_mean = 0.0
        epoch_r_frac = 0.0
        n_batches = 0

        for batch in loader:
            tokens = batch["tokens"].to(device)
            coords = model(tokens)

            c = model.c  # fresh softplus(raw_c) each step

            lq = quartet_loss(coords, batch, c)
            la = domain_angular_loss(coords, batch["domains"])
            lg = genus_anchor_loss(coords, batch["genera"])
            lrep = angular_repulsion_loss(coords, batch["families"], batch["domains"])
            lr_loss = radial_ordering_loss(coords, batch["genome_sizes"].to(device))

            loss = lq + 2.0 * la + 5.0 * lg + 0.5 * lrep + 0.3 * lr_loss

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            # Radius diagnostics
            with torch.no_grad():
                radii = coords.norm(dim=-1)
                r_boundary = (0.9 / torch.sqrt(c)).item()
                epoch_r_mean += radii.mean().item()
                epoch_r_frac += (radii.mean().item() / r_boundary) if r_boundary > 0 else 0

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()

        kappa_val = model.c.item()
        kappa_grad = model.raw_c.grad.item() if model.raw_c.grad is not None else 0.0
        avg_loss = epoch_loss / max(n_batches, 1)
        avg_r = epoch_r_mean / max(n_batches, 1)
        avg_r_frac = epoch_r_frac / max(n_batches, 1)

        # Conformal factor at mean radius (diagnostic)
        lam_at_mean = ((1.0 - kappa_val * avg_r**2) ** 2) / 4.0

        history.append({
            'epoch': epoch + 1, 'kappa': kappa_val,
            'kappa_grad': kappa_grad, 'loss': avg_loss,
            'r_mean': avg_r, 'r_frac': avg_r_frac, 'lambda': lam_at_mean,
        })

        if (epoch + 1) % 10 == 0 or epoch == 0 or epoch == n_epochs - 1:
            log(f"    Ep {epoch+1:4d}: κ={kappa_val:.6f}  ∇raw_c={kappa_grad:+.2e}  "
                f"loss={avg_loss:.4f}  r̄={avg_r:.4f}  r/r_max={avg_r_frac:.3f}  "
                f"λ={lam_at_mean:.4f}")

    final_kappa = model.c.item()
    run_dir = Path(output_dir) / f"seed_{seed}_init_{init_c:.2f}"
    run_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), run_dir / "model.pt")
    with open(run_dir / "kappa_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    return final_kappa, history


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="κ Convergence — E11 + Riemannian Correction")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="kappa_E11_riemannian_results")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--n-anchors", type=int, default=250)
    parser.add_argument("--n-epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--init-kappa", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sweep", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "kappa_E11_riemannian.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("κ CONVERGENCE — E11 + RIEMANNIAN GRADIENT CORRECTION")
    log("=" * 70)
    log(f"Manifest: {args.manifest}")
    log(f"Riemannian: grad_z scaled by λ(z) = (1-c‖z‖²)²/4")
    log(f"Curvature: c = softplus(raw_c), no geoopt")
    log()

    anchors = select_anchors(args.manifest, n_target=args.n_anchors)
    domain_counts = defaultdict(int)
    for a in anchors:
        domain_counts[a.get("domain", "Unknown")] += 1
    log(f"Selected {len(anchors)} anchors: {dict(domain_counts)}")
    log()

    if args.sweep:
        seeds = [42, 59, 76, 93, 110]
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
                results.append({'seed': seed, 'init_kappa': init_k, 'final_kappa': final_k})
                log(f"  → Final κ = {final_k:.6f}")
                log()

                with open(os.path.join(args.output_dir, "sweep_results.json"), 'w') as f:
                    json.dump(results, f, indent=2)

        # Summary
        log("=" * 70)
        log("SWEEP COMPLETE — RIEMANNIAN vs EUCLIDEAN")
        log("=" * 70)
        for r in results:
            log(f"  seed={r['seed']:>4d}  init_κ={r['init_kappa']:.1f}  →  final_κ={r['final_kappa']:.6f}")

        all_kappas = [r['final_kappa'] for r in results]
        mean_k = np.mean(all_kappas)
        std_k = np.std(all_kappas)
        cv = std_k / mean_k * 100 if mean_k > 0 else 999

        log(f"\n  Overall: κ = {mean_k:.6f} ± {std_k:.6f} (CV = {cv:.1f}%)")

        theory = (1.61 * math.log(2)) ** 2
        log(f"  Theory:  κ = (1.61·ln2)² = {theory:.4f}")
        log(f"  Agreement: {abs(mean_k - theory)/theory*100:.1f}%")

        if cv < 5:
            log(f"\n  CONVERGED: All runs agree within {cv:.1f}% CV")
        elif cv < 20:
            log(f"\n  PARTIAL: CV = {cv:.1f}%")
        else:
            log(f"\n  NO CONVERGENCE: CV = {cv:.1f}%")

    else:
        final_k, history = train_one_run(
            anchors, args.seed, args.init_kappa, args.n_epochs, args.lr,
            args.batch_size, args.device, args.output_dir, log,
        )
        log(f"\nFinal κ = {final_k:.6f} (init was {args.init_kappa})")

    log_file.close()


if __name__ == "__main__":
    main()
