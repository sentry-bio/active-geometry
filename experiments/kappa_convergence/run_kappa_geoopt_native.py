#!/usr/bin/env python3
"""
κ Measurement — geoopt-native Direct Manifold Optimization
============================================================

The simplest, most principled κ measurement possible:
  - Each organism is a ManifoldParameter on a PoincareBall
  - Learnable c parameter managed by geoopt (softplus reparameterization)
  - RiemannianAdam for positions (parallel transport, retraction)
  - geoopt.PoincareBall.dist() for all distance computations
  - Same 5 geometric losses from E11

This mirrors how BiosphereCodec works — the same geoopt integration
that produced κ ≈ 1.25 and validated against RNA viruses and proteins.

No manual distance functions. No separate c management. No fighting
the package. Just geoopt doing what it was designed to do.

Key insight: geoopt's in-place softplus reparameterization of c is NOT
corruption — it's a valid parameterization that ensures c > 0 and
provides clean gradient flow during training. The BiosphereCodec's
success proves this works.

Usage:
  python run_kappa_geoopt_native.py --manifest /path/to/manifest.csv --sweep
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
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


# ── Anchor Selection ─────────────────────────────────────────────────────

def select_anchors(manifest_path: str, n_target: int = 500, seed: int = 42) -> List[Dict]:
    """Select taxonomically diverse organisms. No tokenized data needed."""
    rng = random.Random(seed)
    all_rows = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
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

RANK_KEYS = ["genus", "family", "order", "class", "phylum", "domain"]

def taxonomy_distance(a: Dict, b: Dict) -> int:
    for level, key in enumerate(RANK_KEYS):
        va, vb = a.get(key, "Unknown"), b.get(key, "Unknown")
        if va == vb and va not in ("Unknown", ""):
            return level
    return 6


# ── Loss Functions (all using geoopt ball.dist) ──────────────────────────

def quartet_loss(positions: torch.Tensor, ball: geoopt.PoincareBall,
                 organisms: List[Dict], n_quartets: int = 200,
                 margin: float = 0.5) -> torch.Tensor:
    B = positions.shape[0]
    if B < 4:
        return torch.tensor(0.0, device=positions.device)

    indices = list(range(B))
    loss = torch.tensor(0.0, device=positions.device)
    count = 0
    for _ in range(min(n_quartets, B * 3)):
        q = random.sample(indices, 4)
        splits = [(q[0],q[1],q[2],q[3]), (q[0],q[2],q[1],q[3]), (q[0],q[3],q[1],q[2])]
        best = min(splits, key=lambda s:
            taxonomy_distance(organisms[s[0]], organisms[s[1]]) +
            taxonomy_distance(organisms[s[2]], organisms[s[3]]))
        i, j, k, l = best
        d_ij = ball.dist(positions[i], positions[j])
        d_kl = ball.dist(positions[k], positions[l])
        d_ik = ball.dist(positions[i], positions[k])
        d_jl = ball.dist(positions[j], positions[l])
        loss = loss + F.relu((d_ij + d_kl) - (d_ik + d_jl) + margin)
        count += 1
    return loss / max(count, 1)


def domain_angular_loss(positions: torch.Tensor, organisms: List[Dict]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=positions.device)
    count = 0
    for domain, angle_deg in DOMAIN_ANGLES.items():
        rad = math.radians(angle_deg)
        target = torch.tensor([math.cos(rad), math.sin(rad)], device=positions.device)
        for i, org in enumerate(organisms):
            if org.get("domain") == domain:
                norm = positions[i] / positions[i].norm().clamp(min=1e-8)
                loss = loss + (1.0 - (norm * target).sum())
                count += 1
    return loss / max(count, 1)


def genus_anchor_loss(positions: torch.Tensor, organisms: List[Dict]) -> torch.Tensor:
    loss = torch.tensor(0.0, device=positions.device)
    count = 0
    for i, org in enumerate(organisms):
        genus = org.get("genus", "Unknown")
        if genus in GENUS_ANCHORS:
            rad = math.radians(GENUS_ANCHORS[genus])
            target = torch.tensor([math.cos(rad), math.sin(rad)], device=positions.device)
            norm = positions[i] / positions[i].norm().clamp(min=1e-8)
            loss = loss + (1.0 - (norm * target).sum())
            count += 1
    return loss / max(count, 1)


def angular_repulsion_loss(positions: torch.Tensor, organisms: List[Dict],
                            n_pairs: int = 300) -> torch.Tensor:
    B = positions.shape[0]
    if B < 2:
        return torch.tensor(0.0, device=positions.device)
    norms = positions.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    dirs = positions / norms
    loss = torch.tensor(0.0, device=positions.device)
    count = 0
    indices = list(range(B))
    for _ in range(n_pairs):
        i, j = random.sample(indices, 2)
        fi = organisms[i].get("family", "Unknown")
        fj = organisms[j].get("family", "Unknown")
        if fi == fj and fi not in ("Unknown", ""):
            continue
        cos_sim = (dirs[i] * dirs[j]).sum()
        di, dj = organisms[i].get("domain"), organisms[j].get("domain")
        margin = 0.5 if di != dj else 0.966
        loss = loss + F.relu(cos_sim - margin)
        count += 1
    return loss / max(count, 1)


def radial_ordering_loss(positions: torch.Tensor, genome_sizes: torch.Tensor,
                          n_pairs: int = 300, margin: float = 0.02) -> torch.Tensor:
    r = positions.norm(dim=-1)
    valid = (genome_sizes > 0).nonzero(as_tuple=True)[0].tolist()
    if len(valid) < 2:
        return torch.tensor(0.0, device=positions.device)
    loss = torch.tensor(0.0, device=positions.device)
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


# ── Training ─────────────────────────────────────────────────────────────

def train_one_run(
    organisms: List[Dict], seed: int, init_c: float, n_epochs: int,
    lr_positions: float, lr_c: float, device: str, output_dir: str, log,
) -> Tuple[float, List[Dict]]:
    """Direct manifold optimization with geoopt-native geometry."""

    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    N = len(organisms)

    genome_sizes = torch.tensor(
        [float(o.get("genome_size", 0) or 0) for o in organisms],
        dtype=torch.float32, device=device,
    )

    # ── geoopt-native setup ──
    # Learnable curvature as nn.Parameter — geoopt handles reparameterization
    c_param = nn.Parameter(torch.tensor(float(init_c)))
    ball = geoopt.PoincareBall(c=c_param)
    # After init: c_param holds inverse_softplus(init_c)
    # ball.c returns softplus(c_param) = init_c

    log(f"  geoopt reparameterization: raw_c={c_param.item():.6f} → "
        f"ball.c={ball.c.item():.6f}")

    # Organism positions as ManifoldParameter
    init_pos = torch.randn(N, LATENT_DIM, device=device) * 0.1
    positions = geoopt.ManifoldParameter(
        ball.projx(init_pos), manifold=ball,
    )

    # ── Optimizer ──
    # RiemannianAdam for positions (parallel transport + retraction)
    # Regular Adam group for c_param (geoopt handles this correctly —
    # c_param is NOT a ManifoldParameter, so it gets Euclidean updates)
    optimizer = geoopt.optim.RiemannianAdam([
        {'params': [positions], 'lr': lr_positions},
        {'params': [c_param], 'lr': lr_c},
    ], stabilize=10)

    n_params = N * LATENT_DIM + 1
    log(f"  Seed {seed}, init_κ={init_c:.2f}: {N} organisms × {LATENT_DIM}D = "
        f"{n_params:,} params")

    history = []

    for epoch in range(n_epochs):
        # ── Compute distances using geoopt's ball.dist ──
        # ball.c flows through softplus → into dist → into loss → gradients
        c_val = ball.c.item()  # diagnostic only

        lq = quartet_loss(positions, ball, organisms, n_quartets=min(500, N * 2))
        la = domain_angular_loss(positions, organisms)
        lg = genus_anchor_loss(positions, organisms)
        lrep = angular_repulsion_loss(positions, organisms, n_pairs=min(500, N * 2))
        lr_loss = radial_ordering_loss(positions, genome_sizes, n_pairs=min(500, N * 2))

        loss = lq + 2.0 * la + 5.0 * lg + 0.5 * lrep + 0.3 * lr_loss

        optimizer.zero_grad()
        loss.backward()

        # Clip gradients
        torch.nn.utils.clip_grad_norm_([positions, c_param], 1.0)

        optimizer.step()

        # ── Diagnostics ──
        with torch.no_grad():
            kappa = ball.c.item()
            raw_c = c_param.item()
            c_grad = c_param.grad.item() if c_param.grad is not None else 0.0
            radii = positions.norm(dim=-1)
            r_mean = radii.mean().item()
            r_max_obs = radii.max().item()
            r_boundary = (1.0 / math.sqrt(kappa)) if kappa > 1e-7 else 999.0
            r_frac = r_mean / r_boundary if r_boundary > 0 else 0

            # Conformal factor at mean radius
            lam = ((1.0 - kappa * r_mean**2) ** 2) / 4.0

            pos_grad_norm = positions.grad.norm().item() if positions.grad is not None else 0.0

        record = {
            'epoch': epoch + 1, 'kappa': kappa, 'raw_c': raw_c,
            'kappa_grad': c_grad, 'loss': loss.item(),
            'r_mean': r_mean, 'r_max': r_max_obs, 'r_frac': r_frac,
            'lambda': lam, 'pos_grad': pos_grad_norm,
            'lq': lq.item(), 'la': la.item(), 'lg': lg.item(),
            'lrep': lrep.item(), 'lr': lr_loss.item(),
        }
        history.append(record)

        if (epoch + 1) % 10 == 0 or epoch == 0 or epoch == n_epochs - 1:
            log(f"    Ep {epoch+1:4d}: κ={kappa:.6f}  raw_c={raw_c:.4f}  "
                f"∇c={c_grad:+.2e}  loss={loss.item():.4f}  "
                f"r̄={r_mean:.4f}  r/R={r_frac:.3f}  λ={lam:.4f}  "
                f"‖∇pos‖={pos_grad_norm:.4f}")

    final_kappa = ball.c.item()

    # Save
    run_dir = Path(output_dir) / f"seed_{seed}_init_{init_c:.2f}"
    run_dir.mkdir(parents=True, exist_ok=True)
    torch.save({
        'positions': positions.data.cpu(),
        'c_param': c_param.data.cpu(),
        'c_final': final_kappa,
        'ball_c': ball.c.item(),
        'organisms': [{'domain': o.get('domain'), 'phylum': o.get('phylum'),
                       'family': o.get('family'), 'genus': o.get('genus'),
                       'genome_size': o.get('genome_size')} for o in organisms],
    }, run_dir / "state.pt")

    with open(run_dir / "kappa_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    return final_kappa, history


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="κ Measurement — geoopt-native Direct Manifold Optimization")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="kappa_geoopt_native_results")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--n-anchors", type=int, default=500)
    parser.add_argument("--n-epochs", type=int, default=1000,
                        help="More epochs — let κ fully converge")
    parser.add_argument("--lr-positions", type=float, default=1e-2)
    parser.add_argument("--lr-c", type=float, default=1e-2,
                        help="Higher c lr to escape flat regions faster")
    parser.add_argument("--init-kappa", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sweep", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "kappa_geoopt_native.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("κ MEASUREMENT — GEOOPT-NATIVE DIRECT MANIFOLD OPTIMIZATION")
    log("=" * 70)
    log(f"Manifest: {args.manifest}")
    log(f"Method: geoopt.PoincareBall + geoopt.RiemannianAdam")
    log(f"Curvature: geoopt's native softplus reparameterization")
    log(f"Distances: geoopt.PoincareBall.dist() (native)")
    log(f"No encoder — organisms are free ManifoldParameters")
    log(f"This mirrors BiosphereCodec's geoopt integration")
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
                    args.output_dir, log,
                )
                results.append({
                    'seed': seed, 'init_kappa': init_k, 'final_kappa': final_k,
                })
                log(f"  → Final κ = {final_k:.6f}")
                log()

                with open(os.path.join(args.output_dir, "sweep_results.json"), 'w') as f:
                    json.dump(results, f, indent=2)

        # ── Summary ──
        log("=" * 70)
        log("SWEEP COMPLETE — GEOOPT-NATIVE MANIFOLD")
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

        log(f"\n  Convergence by init_κ:")
        for init_k in init_kappas:
            group = [r['final_kappa'] for r in results if r['init_kappa'] == init_k]
            if group:
                log(f"    init={init_k:.2f}: final={np.mean(group):.6f} ± {np.std(group):.6f}")

        with open(os.path.join(args.output_dir, "summary.json"), 'w') as f:
            json.dump({
                'method': 'geoopt_native_riemannian_adam',
                'n_organisms': len(organisms),
                'n_families': len(family_counts),
                'latent_dim': LATENT_DIM,
                'n_epochs': args.n_epochs,
                'lr_positions': args.lr_positions,
                'lr_c': args.lr_c,
                'kappa_mean': float(mean_k),
                'kappa_std': float(std_k),
                'kappa_cv': float(cv),
                'kappa_theory': float(kappa_theory),
                'agreement_pct': float(agreement),
                'results': results,
            }, f, indent=2)

    else:
        log(f"{'='*60}")
        log(f"  SINGLE RUN: seed={args.seed}, init_κ={args.init_kappa}")
        log(f"{'='*60}")

        final_k, history = train_one_run(
            organisms, args.seed, args.init_kappa, args.n_epochs,
            args.lr_positions, args.lr_c, args.device,
            args.output_dir, log,
        )
        log(f"  → Final κ = {final_k:.6f}")

    log_file.close()


if __name__ == "__main__":
    main()
