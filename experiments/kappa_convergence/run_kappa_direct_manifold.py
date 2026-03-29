#!/usr/bin/env python3
"""
κ Measurement via Direct Manifold Optimization
================================================

Strip away the encoder entirely. Each organism is a FREE point on the
Poincaré ball, optimized directly via geoopt's RiemannianAdam.

The question: what curvature κ does the Poincaré ball need such that
organisms — placed solely by taxonomy-driven geometric losses — achieve
optimal separation?

This measures the CONTAINER's intrinsic curvature, not the encoder's
approximation of it.

Architecture:
  - N organism positions as geoopt.ManifoldParameter on PoincareBall
  - RiemannianAdam for positions (proper parallel transport & retraction)
  - Separate Adam for raw_c (softplus reparameterization, avoids geoopt
    c corruption)
  - Manual distance function with learnable c (∂L/∂raw_c flows correctly)
  - Same 5 geometric losses from E11
  - No encoder, no tokenization — pure manifold geometry

Key design:
  - geoopt.PoincareBall is used ONLY for RiemannianAdam's retraction/transport
  - Distance/loss computations use our own c (learnable, in autograd graph)
  - The ball's c is synced from our learnable c each step (slight lag ≈ 0)
  - This isolates the c corruption bug to PoincareBall.__init__ (which we
    work around by passing float values, not tensors)

Usage:
  python run_kappa_direct_manifold.py --manifest /path/to/manifest.csv --sweep
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

import geoopt

# ── Constants ─────────────────────────────────────────────────────────────

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


# ── Manual Poincaré Distance (learnable c in autograd graph) ──────────────

def poincare_distance(x: torch.Tensor, y: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Geodesic distance with learnable c: d(x,y) = (1/√c)·arccosh(...)"""
    eps = 1e-7
    x_sq = (x * x).sum(-1, keepdim=True)
    y_sq = (y * y).sum(-1, keepdim=True)
    xy_diff_sq = ((x - y) ** 2).sum(-1, keepdim=True)
    num = 2.0 * c * xy_diff_sq
    denom = (1.0 - c * x_sq).clamp(min=eps) * (1.0 - c * y_sq).clamp(min=eps)
    arg = 1.0 + num / denom
    return ((1.0 / torch.sqrt(c + eps)) * torch.acosh(arg.clamp(min=1.0 + eps))).squeeze(-1)


def conformal_factor(z: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """λ(z) = (1 - c‖z‖²)² / 4  — diagnostic only."""
    z_sq = (z * z).sum(-1, keepdim=True)
    return ((1.0 - c * z_sq).clamp(min=1e-7) ** 2) / 4.0


# ── Anchor Selection (no tokenized paths needed!) ────────────────────────

def select_anchors(manifest_path: str, n_target: int = 500, seed: int = 42) -> List[Dict]:
    """Select taxonomically diverse organisms. No tokenized data needed."""
    rng = random.Random(seed)
    all_rows = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            if row.get("domain") in ("Viruses", ""):
                continue
            # We only need taxonomy + genome_size — no tokenized path required
            all_rows.append(row)

    # Group by domain/phylum/family for stratified sampling
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

    # Ensure genus anchors are present
    for genus in GENUS_ANCHORS:
        if not any(a.get("genus") == genus for a in anchors):
            for row in all_rows:
                if row.get("genus") == genus:
                    anchors.append(row)
                    break

    # Trim if too many
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

RANK_KEYS = ["genus", "family", "order", "class", "phylum", "domain"]

def taxonomy_distance(a: Dict, b: Dict) -> int:
    """How many ranks up until a and b share a label?
    0 = same genus, 1 = same family, ..., 6 = no match."""
    for level, key in enumerate(RANK_KEYS):
        va, vb = a.get(key, "Unknown"), b.get(key, "Unknown")
        if va == vb and va not in ("Unknown", ""):
            return level
    return 6


def batch_taxonomy_distance(organisms: List[Dict], i: int, j: int) -> int:
    return taxonomy_distance(organisms[i], organisms[j])


# ── Loss Functions (same 5 from E11, adapted for direct positions) ────────

def quartet_loss(coords: torch.Tensor, organisms: List[Dict], c: torch.Tensor,
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
        best = min(splits, key=lambda s:
            batch_taxonomy_distance(organisms, s[0], s[1]) +
            batch_taxonomy_distance(organisms, s[2], s[3]))
        i, j, k, l = best
        d_ij = poincare_distance(coords[i:i+1], coords[j:j+1], c)
        d_kl = poincare_distance(coords[k:k+1], coords[l:l+1], c)
        d_ik = poincare_distance(coords[i:i+1], coords[k:k+1], c)
        d_jl = poincare_distance(coords[j:j+1], coords[l:l+1], c)
        loss = loss + F.relu((d_ij + d_kl) - (d_ik + d_jl) + margin).squeeze()
        count += 1
    return loss / max(count, 1)


def domain_angular_loss(coords: torch.Tensor, organisms: List[Dict]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for domain, angle_deg in DOMAIN_ANGLES.items():
        rad = math.radians(angle_deg)
        target = torch.tensor([math.cos(rad), math.sin(rad)], device=coords.device)
        for i, org in enumerate(organisms):
            if org.get("domain") == domain:
                norm = coords[i] / coords[i].norm().clamp(min=1e-8)
                loss = loss + (1.0 - (norm * target).sum())
                count += 1
    return loss / max(count, 1)


def genus_anchor_loss(coords: torch.Tensor, organisms: List[Dict]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=coords.device)
    count = 0
    for i, org in enumerate(organisms):
        genus = org.get("genus", "Unknown")
        if genus in GENUS_ANCHORS:
            rad = math.radians(GENUS_ANCHORS[genus])
            target = torch.tensor([math.cos(rad), math.sin(rad)], device=coords.device)
            norm = coords[i] / coords[i].norm().clamp(min=1e-8)
            loss = loss + (1.0 - (norm * target).sum())
            count += 1
    return loss / max(count, 1)


def angular_repulsion_loss(coords: torch.Tensor, organisms: List[Dict],
                            n_pairs: int = 300) -> torch.Tensor:
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
        fi, fj = organisms[i].get("family", "Unknown"), organisms[j].get("family", "Unknown")
        if fi == fj and fi not in ("Unknown", ""):
            continue
        cos_sim = (dirs[i] * dirs[j]).sum()
        di, dj = organisms[i].get("domain"), organisms[j].get("domain")
        margin = 0.5 if di != dj else 0.966
        loss = loss + F.relu(cos_sim - margin)
        count += 1
    return loss / max(count, 1)


def radial_ordering_loss(coords: torch.Tensor, genome_sizes: torch.Tensor,
                          n_pairs: int = 300, margin: float = 0.02) -> torch.Tensor:
    r = coords.norm(dim=-1)
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


# ── Direct Manifold Optimization ─────────────────────────────────────────

def train_one_run(
    organisms: List[Dict], seed: int, init_c: float, n_epochs: int,
    lr_positions: float, lr_c: float, device: str, output_dir: str, log,
    n_organisms: int = None,
) -> Tuple[float, List[Dict]]:
    """Optimize organism positions directly on the Poincaré ball."""

    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    N = len(organisms) if n_organisms is None else min(n_organisms, len(organisms))
    orgs = organisms[:N]

    # Genome sizes for radial ordering
    genome_sizes = torch.tensor(
        [float(o.get("genome_size", 0) or 0) for o in orgs],
        dtype=torch.float32, device=device,
    )

    # ── Curvature: softplus reparameterization (our own, not geoopt's) ──
    if init_c > 0.01:
        raw_init = math.log(math.exp(init_c) - 1.0)
    else:
        raw_init = -4.0
    raw_c = nn.Parameter(torch.tensor(raw_init, device=device))

    def get_c():
        return F.softplus(raw_c)

    # ── Organism positions: ManifoldParameter on PoincareBall ──
    # Initialize with small random positions (near origin)
    init_positions = torch.randn(N, LATENT_DIM, device=device) * 0.1

    # Create PoincareBall with FLOAT value of c (avoids tensor corruption)
    ball = geoopt.PoincareBall(c=init_c)
    positions = geoopt.ManifoldParameter(
        ball.projx(init_positions),
        manifold=ball,
    )

    # ── Optimizers ──
    # RiemannianAdam for positions: proper parallel transport & retraction
    riemannian_opt = geoopt.optim.RiemannianAdam(
        [positions], lr=lr_positions, stabilize=10,
    )
    # Regular Adam for curvature
    c_opt = torch.optim.Adam([raw_c], lr=lr_c)

    n_params = N * LATENT_DIM + 1  # positions + raw_c
    log(f"  Seed {seed}, init_κ={init_c:.2f}: {N} organisms × {LATENT_DIM}D = "
        f"{n_params:,} params, raw_c={raw_c.item():.4f}, c={get_c().item():.6f}")

    history = []

    for epoch in range(n_epochs):
        # ── Sync ball curvature with current c ──
        # geoopt uses this for retraction/transport in RiemannianAdam
        c_val = get_c().item()
        positions.manifold = geoopt.PoincareBall(c=c_val)

        # ── Forward: all organisms at once (no batching needed) ──
        c = get_c()  # fresh softplus(raw_c) in autograd graph
        coords = positions  # these ARE the positions (ManifoldParameter)

        # ── 5-term geometric loss ──
        lq = quartet_loss(coords, orgs, c, n_quartets=min(500, N * 2))
        la = domain_angular_loss(coords, orgs)
        lg = genus_anchor_loss(coords, orgs)
        lrep = angular_repulsion_loss(coords, orgs, n_pairs=min(500, N * 2))
        lr_loss = radial_ordering_loss(coords, genome_sizes, n_pairs=min(500, N * 2))

        loss = lq + 2.0 * la + 5.0 * lg + 0.5 * lrep + 0.3 * lr_loss

        # ── Backward ──
        riemannian_opt.zero_grad()
        c_opt.zero_grad()
        loss.backward()

        # Clip gradients
        torch.nn.utils.clip_grad_norm_([positions], 1.0)
        torch.nn.utils.clip_grad_norm_([raw_c], 1.0)

        # ── Step: RiemannianAdam for positions, Adam for c ──
        riemannian_opt.step()
        c_opt.step()

        # ── Diagnostics ──
        with torch.no_grad():
            kappa = get_c().item()
            c_grad = raw_c.grad.item() if raw_c.grad is not None else 0.0
            radii = positions.norm(dim=-1)
            r_mean = radii.mean().item()
            r_max = radii.max().item()
            r_boundary = 0.9 / math.sqrt(kappa) if kappa > 1e-7 else 999.0
            r_frac = r_mean / r_boundary if r_boundary > 0 else 0
            lam = conformal_factor(positions, get_c()).mean().item()

            # Position gradient magnitude (Riemannian)
            pos_grad_norm = positions.grad.norm().item() if positions.grad is not None else 0.0

        record = {
            'epoch': epoch + 1, 'kappa': kappa,
            'kappa_grad': c_grad, 'loss': loss.item(),
            'r_mean': r_mean, 'r_max': r_max, 'r_frac': r_frac,
            'lambda': lam, 'pos_grad': pos_grad_norm,
            'lq': lq.item(), 'la': la.item(), 'lg': lg.item(),
            'lrep': lrep.item(), 'lr': lr_loss.item(),
        }
        history.append(record)

        if (epoch + 1) % 10 == 0 or epoch == 0 or epoch == n_epochs - 1:
            log(f"    Ep {epoch+1:4d}: κ={kappa:.6f}  ∇raw_c={c_grad:+.2e}  "
                f"loss={loss.item():.4f}  r̄={r_mean:.4f}  r/r_max={r_frac:.3f}  "
                f"λ={lam:.4f}  ‖∇pos‖={pos_grad_norm:.4f}")

    final_kappa = get_c().item()
    run_dir = Path(output_dir) / f"seed_{seed}_init_{init_c:.2f}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Save positions and curvature
    torch.save({
        'positions': positions.data.cpu(),
        'raw_c': raw_c.data.cpu(),
        'c_final': final_kappa,
        'organisms': [{'domain': o.get('domain'), 'phylum': o.get('phylum'),
                       'family': o.get('family'), 'genus': o.get('genus'),
                       'genome_size': o.get('genome_size')} for o in orgs],
    }, run_dir / "state.pt")

    with open(run_dir / "kappa_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    return final_kappa, history


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="κ Measurement — Direct Manifold Optimization (geoopt RiemannianAdam)")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="kappa_direct_manifold_results")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--n-organisms", type=int, default=None,
                        help="Max organisms (default: all selected)")
    parser.add_argument("--n-anchors", type=int, default=500,
                        help="Target anchor count for selection")
    parser.add_argument("--n-epochs", type=int, default=500)
    parser.add_argument("--lr-positions", type=float, default=1e-2,
                        help="Learning rate for organism positions (RiemannianAdam)")
    parser.add_argument("--lr-c", type=float, default=3e-3,
                        help="Learning rate for curvature (Adam)")
    parser.add_argument("--init-kappa", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sweep", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "kappa_direct_manifold.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("κ MEASUREMENT — DIRECT MANIFOLD OPTIMIZATION")
    log("=" * 70)
    log(f"Manifest: {args.manifest}")
    log(f"Method: geoopt.RiemannianAdam for positions on PoincareBall")
    log(f"Curvature: c = softplus(raw_c), managed outside geoopt")
    log(f"No encoder — organisms are free ManifoldParameters")
    log()

    organisms = select_anchors(args.manifest, n_target=args.n_anchors)
    domain_counts = defaultdict(int)
    family_counts = set()
    for o in organisms:
        domain_counts[o.get("domain", "Unknown")] += 1
        family_counts.add(o.get("family", "Unknown"))
    log(f"Selected {len(organisms)} organisms: {dict(domain_counts)}")
    log(f"  {len(family_counts)} unique families")
    log()

    if args.sweep:
        seeds = [42, 59, 76, 93, 110]
        init_kappas = [0.1, 0.5, 1.0, 1.25, 2.0, 5.0]
        results = []

        for seed in seeds:
            for init_k in init_kappas:
                log(f"{'='*60}")
                log(f"  SWEEP: seed={seed}, init_κ={init_k}")
                log(f"{'='*60}")

                final_k, history = train_one_run(
                    organisms, seed, init_k, args.n_epochs,
                    args.lr_positions, args.lr_c, args.device,
                    args.output_dir, log, args.n_organisms,
                )
                results.append({
                    'seed': seed, 'init_kappa': init_k, 'final_kappa': final_k,
                })
                log(f"  → Final κ = {final_k:.6f}")
                log()

                # Save incrementally
                with open(os.path.join(args.output_dir, "sweep_results.json"), 'w') as f:
                    json.dump(results, f, indent=2)

        # ── Summary ──
        log("=" * 70)
        log("SWEEP COMPLETE — DIRECT MANIFOLD OPTIMIZATION")
        log("=" * 70)
        log(f"  {'seed':>6s}  {'init_κ':>8s}  {'final_κ':>10s}")
        log(f"  {'─'*6}  {'─'*8}  {'─'*10}")
        for r in results:
            log(f"  {r['seed']:>6d}  {r['init_kappa']:>8.2f}  {r['final_kappa']:>10.6f}")

        all_kappas = [r['final_kappa'] for r in results]
        mean_k = np.mean(all_kappas)
        std_k = np.std(all_kappas)
        cv = std_k / mean_k * 100 if mean_k > 0 else 999

        kappa_theory = (1.6 * math.log(2)) ** 2
        agreement = abs(mean_k - kappa_theory) / kappa_theory * 100

        log(f"\n  κ = {mean_k:.6f} ± {std_k:.6f} (CV = {cv:.1f}%)")
        log(f"  Theory: κ = (1.6·ln2)² = {kappa_theory:.6f}")
        log(f"  Agreement: {agreement:.1f}%")

        # Group by init_kappa to check convergence
        log(f"\n  Convergence by init_κ:")
        for init_k in init_kappas:
            group = [r['final_kappa'] for r in results if r['init_kappa'] == init_k]
            if group:
                log(f"    init={init_k:.2f}: final={np.mean(group):.6f} ± {np.std(group):.6f}")

        with open(os.path.join(args.output_dir, "summary.json"), 'w') as f:
            json.dump({
                'method': 'direct_manifold_geoopt_riemannian_adam',
                'n_organisms': len(organisms),
                'n_families': len(family_counts),
                'latent_dim': LATENT_DIM,
                'kappa_mean': float(mean_k),
                'kappa_std': float(std_k),
                'kappa_cv': float(cv),
                'kappa_theory': float(kappa_theory),
                'agreement_pct': float(agreement),
                'results': results,
            }, f, indent=2)

    else:
        # Single run
        log(f"{'='*60}")
        log(f"  SINGLE RUN: seed={args.seed}, init_κ={args.init_kappa}")
        log(f"{'='*60}")

        final_k, history = train_one_run(
            organisms, args.seed, args.init_kappa, args.n_epochs,
            args.lr_positions, args.lr_c, args.device,
            args.output_dir, log, args.n_organisms,
        )
        log(f"  → Final κ = {final_k:.6f}")

    log_file.close()


if __name__ == "__main__":
    main()
