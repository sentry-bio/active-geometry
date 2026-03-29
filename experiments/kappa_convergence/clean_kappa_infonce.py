#!/usr/bin/env python3
"""
clean_kappa_infonce.py — Clean κ measurement via InfoNCE

The simplest possible experiment to measure the optimal Poincaré ball
curvature for genomic embeddings. No geoopt, no prior, no complexity.

Architecture:
  Embedding(4096, d_model) → mean_pool → Linear(d_model, 2) → exp_map → H²

Loss:
  InfoNCE with Poincaré distances: same-genus = positive, different = negative
  κ is load-bearing because distance ratios depend on curvature.

Protocol:
  5 curvature initializations × 3 seeds = 15 runs
  Track c = softplus(raw_c) every epoch
  If all converge to same c → that's the measured κ
"""

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader


# ─── Manual Poincaré geometry (no geoopt) ─────────────────────────────

def poincare_distance(u, v, c, eps=1e-7):
    """Poincaré ball distance: d(u,v) = (1/√c) · acosh(1 + 2c‖u-v‖²/((1-c‖u‖²)(1-c‖v‖²)))"""
    diff_sq = torch.sum((u - v) ** 2, dim=-1)
    u_sq = torch.sum(u ** 2, dim=-1)
    v_sq = torch.sum(v ** 2, dim=-1)
    denom = (1.0 - c * u_sq).clamp(min=eps) * (1.0 - c * v_sq).clamp(min=eps)
    arg = 1.0 + 2.0 * c * diff_sq / denom
    return (1.0 / torch.sqrt(c + eps)) * torch.acosh(arg.clamp(min=1.0 + eps))


def exp_map_origin(v, c, eps=1e-7):
    """Exponential map from origin in Poincaré ball with curvature c."""
    sqrt_c = torch.sqrt(c + eps)
    v_norm = torch.norm(v, dim=-1, keepdim=True).clamp(min=eps)
    return torch.tanh(sqrt_c * v_norm / 2.0) * v / (sqrt_c * v_norm)


def project_to_ball(x, c, max_norm_frac=0.95):
    """Project points to stay inside the Poincaré ball of radius 1/√c."""
    r_max = max_norm_frac / torch.sqrt(c + 1e-7)
    norms = torch.norm(x, dim=-1, keepdim=True)
    scale = torch.clamp(r_max / (norms + 1e-7), max=1.0)
    return x * scale


# ─── Minimal encoder ──────────────────────────────────────────────────

class MinimalPoincareEncoder(nn.Module):
    """
    Tiny encoder: Embed tokens → mean pool → project to H².
    ~100K parameters. κ is the only thing that matters.
    """
    def __init__(self, vocab_size=4096, d_model=64, d_hyp=2, c_init=1.0):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, d_model, padding_idx=0)
        self.proj = nn.Linear(d_model, d_hyp)
        # Learnable curvature via softplus reparameterization
        # c = softplus(raw_c), so raw_c = softplus_inverse(c_init)
        raw_c_init = math.log(math.exp(c_init) - 1) if c_init > 0.01 else -4.0
        self.raw_c = nn.Parameter(torch.tensor(raw_c_init))

    @property
    def c(self):
        """Actual curvature used in all computations."""
        return F.softplus(self.raw_c)

    def forward(self, tokens):
        """tokens: [B, L] → embeddings in H²: [B, d_hyp]"""
        x = self.embed(tokens)            # [B, L, d_model]
        x = x.mean(dim=1)                 # [B, d_model]
        x = self.proj(x)                  # [B, d_hyp]
        c = self.c
        x = exp_map_origin(x, c)          # [B, d_hyp] on Poincaré ball
        x = project_to_ball(x, c)         # ensure inside ball
        return x


# ─── InfoNCE loss with Poincaré distances ─────────────────────────────

def infonce_poincare(embeddings, labels, c, temperature=0.1):
    """
    InfoNCE using Poincaré distances.

    For each anchor, positives = same genus, negatives = different genus.
    loss = -log(exp(-d(a,p)/τ) / Σ_j exp(-d(a,j)/τ))

    κ is load-bearing: changing c changes all distance ratios.
    """
    B = embeddings.shape[0]

    # Pairwise distances [B, B]
    u = embeddings.unsqueeze(1).expand(B, B, -1)  # [B, B, d]
    v = embeddings.unsqueeze(0).expand(B, B, -1)  # [B, B, d]
    D = poincare_distance(u, v, c)                 # [B, B]

    # Similarity = negative distance / temperature
    sim = -D / temperature                         # [B, B]

    # Mask: same genus = positive
    label_eq = labels.unsqueeze(0) == labels.unsqueeze(1)  # [B, B]
    # Exclude self-pairs
    eye_mask = ~torch.eye(B, dtype=torch.bool, device=embeddings.device)
    pos_mask = label_eq & eye_mask
    neg_mask = ~label_eq & eye_mask

    # Skip anchors with no positives or no negatives
    has_pos = pos_mask.any(dim=1)
    has_neg = neg_mask.any(dim=1)
    valid = has_pos & has_neg

    if valid.sum() == 0:
        return torch.tensor(0.0, device=embeddings.device, requires_grad=True)

    # For each valid anchor: log_sum_exp over all non-self, minus log_sum_exp over positives
    # Standard supervised InfoNCE
    sim_valid = sim[valid]                         # [V, B]
    pos_valid = pos_mask[valid]                    # [V, B]
    all_valid = eye_mask[valid]                    # [V, B]

    # Denominator: sum over all non-self
    sim_all = sim_valid.clone()
    sim_all[~all_valid] = -1e9
    log_denom = torch.logsumexp(sim_all, dim=1)   # [V]

    # Numerator: mean over positives (for each anchor, average positive log-prob)
    loss = 0.0
    for i in range(sim_valid.shape[0]):
        pos_sims = sim_valid[i][pos_valid[i]]     # positives for anchor i
        log_num = torch.logsumexp(pos_sims, dim=0) - math.log(pos_sims.shape[0])
        loss += -(log_num - log_denom[i])

    return loss / sim_valid.shape[0]


# ─── Dataset ──────────────────────────────────────────────────────────

class GenomeDataset(Dataset):
    def __init__(self, manifest_path, max_genomes=None, min_genus_count=5):
        import pandas as pd
        df = pd.read_csv(manifest_path)
        df = df[df['tokenized_path'].notna()]
        df = df[df['genus'].notna()]
        df = df[df['tokenized_path'].apply(lambda p: os.path.exists(str(p)))]

        # Filter genera with at least min_genus_count members (need positives)
        gc = df['genus'].value_counts()
        valid_genera = gc[gc >= min_genus_count].index
        df = df[df['genus'].isin(valid_genera)]

        if max_genomes and len(df) > max_genomes:
            df = df.sample(max_genomes, random_state=42)

        self.paths = df['tokenized_path'].tolist()
        # Encode genera as integers
        genera = df['genus'].tolist()
        genus_to_id = {g: i for i, g in enumerate(sorted(set(genera)))}
        self.labels = [genus_to_id[g] for g in genera]
        self.n_genera = len(genus_to_id)

        print(f"Dataset: {len(self.paths)} genomes, {self.n_genera} genera")

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        tokens = np.load(self.paths[idx]).astype(np.int64)
        return torch.from_numpy(tokens), self.labels[idx]


# ─── Training ─────────────────────────────────────────────────────────

def train_single_run(config):
    """Train one run with given c_init and seed. Returns curvature trajectory."""
    torch.manual_seed(config['seed'])
    np.random.seed(config['seed'])

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    dataset = GenomeDataset(
        config['manifest_path'],
        max_genomes=config.get('max_genomes'),
        min_genus_count=config.get('min_genus_count', 5),
    )

    loader = DataLoader(
        dataset, batch_size=config['batch_size'], shuffle=True,
        num_workers=0, drop_last=True,
    )

    model = MinimalPoincareEncoder(
        vocab_size=config.get('vocab_size', 4096),
        d_model=config.get('d_model', 64),
        d_hyp=config.get('d_hyp', 2),
        c_init=config['c_init'],
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"\nRun: c_init={config['c_init']:.2f}, seed={config['seed']}, "
          f"params={n_params:,}, device={device}")
    print(f"  c(0) = {model.c.item():.6f} (raw_c = {model.raw_c.item():.4f})")

    # Separate optimizer groups: encoder params + curvature
    encoder_params = [p for n, p in model.named_parameters() if n != 'raw_c']
    optimizer = torch.optim.Adam([
        {'params': encoder_params, 'lr': config.get('lr', 1e-3)},
        {'params': [model.raw_c], 'lr': config.get('c_lr', 1e-2), 'weight_decay': 0},
    ])

    history = []
    t0 = time.time()

    for epoch in range(config['epochs']):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for tokens, labels in loader:
            tokens = tokens.to(device)
            labels = torch.tensor(labels, device=device)

            optimizer.zero_grad()
            emb = model(tokens)
            loss = infonce_poincare(
                emb, labels, model.c,
                temperature=config.get('temperature', 0.1),
            )

            if loss.item() > 0:
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        c_val = model.c.item()
        raw_c_val = model.raw_c.item()

        # Check gradient on curvature
        c_grad = model.raw_c.grad.item() if model.raw_c.grad is not None else 0.0

        record = {
            'epoch': epoch,
            'loss': avg_loss,
            'c': c_val,
            'raw_c': raw_c_val,
            'c_grad': c_grad,
            'time': time.time() - t0,
        }
        history.append(record)

        if epoch % config.get('log_every', 5) == 0 or epoch == config['epochs'] - 1:
            print(f"  Epoch {epoch:3d} | loss={avg_loss:.4f} | "
                  f"c={c_val:.6f} | ∇c={c_grad:.6f} | "
                  f"t={time.time()-t0:.0f}s")

    return {
        'c_init': config['c_init'],
        'seed': config['seed'],
        'c_final': model.c.item(),
        'loss_final': history[-1]['loss'],
        'history': history,
        'n_params': n_params,
        'n_genomes': len(dataset),
        'n_genera': dataset.n_genera,
    }


# ─── Main sweep ───────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Clean κ measurement via InfoNCE")
    parser.add_argument('--manifest', type=str, required=True,
                        help='Path to manifest CSV with tokenized_path and genus columns')
    parser.add_argument('--output', type=str, default='clean_kappa_results',
                        help='Output directory')
    parser.add_argument('--max-genomes', type=int, default=None,
                        help='Max genomes to use (None = all)')
    parser.add_argument('--min-genus-count', type=int, default=5,
                        help='Min genomes per genus to include')
    parser.add_argument('--epochs', type=int, default=50,
                        help='Training epochs per run')
    parser.add_argument('--batch-size', type=int, default=128,
                        help='Batch size')
    parser.add_argument('--d-model', type=int, default=64,
                        help='Embedding dimension')
    parser.add_argument('--d-hyp', type=int, default=2,
                        help='Poincaré ball dimension')
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--c-lr', type=float, default=1e-2,
                        help='Learning rate for curvature parameter')
    parser.add_argument('--temperature', type=float, default=0.1,
                        help='InfoNCE temperature')
    parser.add_argument('--c-inits', type=str, default='0.1,0.5,1.0,2.0,5.0',
                        help='Comma-separated c_init values to sweep')
    parser.add_argument('--seeds', type=str, default='42,137,2024',
                        help='Comma-separated random seeds')
    parser.add_argument('--log-every', type=int, default=5)
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: 500 genomes, 10 epochs')
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    c_inits = [float(x) for x in args.c_inits.split(',')]
    seeds = [int(x) for x in args.seeds.split(',')]

    if args.quick:
        args.max_genomes = 500
        args.epochs = 10
        args.log_every = 1

    base_config = {
        'manifest_path': args.manifest,
        'max_genomes': args.max_genomes,
        'min_genus_count': args.min_genus_count,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'd_model': args.d_model,
        'd_hyp': args.d_hyp,
        'lr': args.lr,
        'c_lr': args.c_lr,
        'temperature': args.temperature,
        'log_every': args.log_every,
    }

    print("=" * 60)
    print("CLEAN κ MEASUREMENT — InfoNCE + Manual Poincaré Distance")
    print("=" * 60)
    print(f"No geoopt. No prior. No complexity.")
    print(f"c_inits: {c_inits}")
    print(f"seeds:   {seeds}")
    print(f"Total runs: {len(c_inits) * len(seeds)}")
    print()

    all_results = []

    for c_init in c_inits:
        for seed in seeds:
            config = {**base_config, 'c_init': c_init, 'seed': seed}
            result = train_single_run(config)
            all_results.append(result)

            # Save incrementally
            with open(outdir / 'results.json', 'w') as f:
                json.dump({
                    'config': base_config,
                    'c_inits': c_inits,
                    'seeds': seeds,
                    'runs': [{k: v for k, v in r.items() if k != 'history'}
                             for r in all_results],
                }, f, indent=2)

            # Save full history per run
            run_name = f"c{c_init:.2f}_s{seed}"
            with open(outdir / f'{run_name}_history.json', 'w') as f:
                json.dump(result['history'], f, indent=2)

    # ─── Analysis ───────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"{'c_init':>8s} {'seed':>6s} {'c_final':>10s} {'loss':>8s}")
    print(f"{'─'*8} {'─'*6} {'─'*10} {'─'*8}")

    c_finals = []
    for r in all_results:
        print(f"{r['c_init']:>8.2f} {r['seed']:>6d} {r['c_final']:>10.6f} {r['loss_final']:>8.4f}")
        c_finals.append(r['c_final'])

    c_finals = np.array(c_finals)
    print(f"\n{'─'*40}")
    print(f"Mean c_final:   {c_finals.mean():.6f}")
    print(f"Std c_final:    {c_finals.std():.6f}")
    print(f"CV:             {c_finals.std()/c_finals.mean()*100:.2f}%")
    print(f"Min:            {c_finals.min():.6f}")
    print(f"Max:            {c_finals.max():.6f}")

    # Check convergence: do all runs agree?
    cv = c_finals.std() / c_finals.mean() * 100
    if cv < 5:
        print(f"\n✓ CONVERGED: All runs agree within {cv:.1f}% CV")
        print(f"  κ = {c_finals.mean():.4f} ± {c_finals.std():.4f}")

        # Compare with theory
        h = 1.6  # bits/nt
        kappa_theory = (h * math.log(2)) ** 2
        agreement = abs(c_finals.mean() - kappa_theory) / kappa_theory * 100
        print(f"\n  Theory prediction: κ = (h·ln2)² = ({h}×{math.log(2):.4f})² = {kappa_theory:.4f}")
        print(f"  Agreement: {agreement:.1f}%")
    elif cv < 20:
        print(f"\n⚠ PARTIAL CONVERGENCE: CV = {cv:.1f}%")
        # Group by c_init to check if init matters
        for ci in sorted(set(r['c_init'] for r in all_results)):
            subset = [r['c_final'] for r in all_results if r['c_init'] == ci]
            print(f"  c_init={ci:.2f}: c_final = {np.mean(subset):.6f} ± {np.std(subset):.6f}")
    else:
        print(f"\n✗ NO CONVERGENCE: CV = {cv:.1f}%")
        print(f"  κ is initialization-dependent — no universal attractor found")

    # Save final summary
    summary = {
        'config': base_config,
        'c_inits': c_inits,
        'seeds': seeds,
        'c_final_mean': float(c_finals.mean()),
        'c_final_std': float(c_finals.std()),
        'c_final_cv_pct': float(cv),
        'c_finals': c_finals.tolist(),
        'runs': [{k: v for k, v in r.items() if k != 'history'} for r in all_results],
        'theory': {
            'h': 1.6,
            'kappa_predicted': float(kappa_theory) if cv < 20 else None,
            'agreement_pct': float(agreement) if cv < 5 else None,
        },
    }
    with open(outdir / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nAll results saved to {outdir}/")


if __name__ == '__main__':
    main()
