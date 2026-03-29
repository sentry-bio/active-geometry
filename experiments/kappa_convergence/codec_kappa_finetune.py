#!/usr/bin/env python3
"""
codec_kappa_finetune.py — Measure κ using pretrained BiosphereCodec

Loads a pretrained checkpoint, patches PoincareMapping to use manual
Poincaré distance (no geoopt), sets κ learnable, and fine-tunes with
InfoNCE active (tax_ids passed to forward).

The encoder already learned to compress genomes (h is embedded in the
representations). This script asks: what curvature κ best organizes
those compressed representations in H²?

Supports --project-dim to project 128D representations down to lower
dimensions (e.g., 2D). When set, a learnable linear projection is
added after the encoder, and κ is measured in the projected space.
This lets us test whether the κ = 5/4 signal survives projection to
the 2D Poincaré disk H².

No geoopt. No prior. Compression already done by pretraining.
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader


# ─── Manual Poincaré geometry ─────────────────────────────────────────

def manual_poincare_distance(u, v, c, eps=1e-7):
    """Standard Poincaré ball distance with curvature c."""
    diff_sq = torch.sum((u - v) ** 2, dim=-1)
    u_sq = torch.sum(u ** 2, dim=-1)
    v_sq = torch.sum(v ** 2, dim=-1)
    denom = (1.0 - c * u_sq).clamp(min=eps) * (1.0 - c * v_sq).clamp(min=eps)
    arg = 1.0 + 2.0 * c * diff_sq / denom
    return (1.0 / torch.sqrt(c + eps)) * torch.acosh(arg.clamp(min=1.0 + eps))


def manual_projx(x, c, max_norm_frac=0.95):
    """Project to Poincaré ball interior."""
    r_max = max_norm_frac / torch.sqrt(c + 1e-7)
    norms = torch.norm(x, dim=-1, keepdim=True)
    scale = torch.clamp(r_max / (norms + 1e-7), max=1.0)
    return x * scale


# ─── Monkey-patch PoincareMapping ─────────────────────────────────────

def patch_poincare_mapping(hyper_module):
    """
    Replace geoopt-based methods with manual Poincaré math.
    This avoids the in-place softplus corruption entirely.
    """
    # Store original c value (may already be corrupted by geoopt init)
    # Check if geoopt already corrupted it
    original_c = hyper_module.c.item()

    def patched_forward(self, x):
        z_euc = torch.tanh(self.lin(x))
        r_max = 0.9 / torch.sqrt(self.c.clamp(min=1e-7))
        z_euc = z_euc * r_max
        return manual_projx(z_euc, self.c)

    def patched_dist_mat(self, z):
        B = z.shape[0]
        u = z.unsqueeze(1).expand(B, B, -1)
        v = z.unsqueeze(0).expand(B, B, -1)
        return manual_poincare_distance(u, v, self.c)

    import types
    hyper_module.forward = types.MethodType(patched_forward, hyper_module)
    hyper_module.dist_mat = types.MethodType(patched_dist_mat, hyper_module)
    # Kill the geoopt cache to prevent any further geoopt calls
    hyper_module._manifold = None
    hyper_module._c_cached = None
    # Override _man to return None (forces fallback paths if any code calls it)
    hyper_module._man = lambda: None

    return original_c


# ─── Dataset ──────────────────────────────────────────────────────────

class TokenizedGenomeDataset(Dataset):
    """Load tokenized genomes with taxonomy labels for InfoNCE."""

    def __init__(self, manifest_path, max_genomes=None, min_genus_count=2):
        import pandas as pd
        df = pd.read_csv(manifest_path, low_memory=False)
        df = df[df['tokenized_path'].notna() & df['genus'].notna()]
        df = df[df['tokenized_path'].apply(lambda p: os.path.exists(str(p)))]

        # Need at least min_genus_count per genus for InfoNCE positives
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

        print(f"Dataset: {len(self.paths)} genomes, {self.n_genera} genera")

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        tokens = np.load(self.paths[idx]).astype(np.int64)
        return torch.from_numpy(tokens), self.labels[idx]


# ─── Training ─────────────────────────────────────────────────────────

def run_finetune(config):
    """Fine-tune pretrained BiosphereCodec with learnable κ."""
    torch.manual_seed(config['seed'])
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # --- Load BiosphereCodec ---
    sys.path.insert(0, config['codec_dir'])
    from BiosphereCodec import BiosphereCodec

    model = BiosphereCodec(
        vocab=config['vocab_size'],
        d_model=config['d_model'],
        n_layers=config.get('n_layers', 4),
        max_len=config.get('max_len', 8192),
        latent_dim=config['latent_dim'],
    )

    # Load pretrained weights
    checkpoint = torch.load(config['checkpoint_path'], map_location='cpu', weights_only=False)
    state = checkpoint.get('model_state_dict', checkpoint)
    missing, unexpected = model.load_state_dict(state, strict=False)
    print(f"Loaded checkpoint: {len(state)} keys, {len(missing)} missing, {len(unexpected)} unexpected")

    # --- Read the actual κ from checkpoint ---
    c_from_checkpoint = model.hyper.c.item()
    print(f"κ from checkpoint (raw parameter): {c_from_checkpoint:.6f}")

    # --- Patch to bypass geoopt ---
    original_c = patch_poincare_mapping(model.hyper)
    print(f"Patched PoincareMapping: manual Poincaré distance (no geoopt)")

    # --- Set κ to desired initial value ---
    c_init = config['c_init']
    model.hyper.c.data.fill_(c_init)
    model.hyper.c.requires_grad = True
    print(f"κ set to {c_init:.4f} (learnable)")

    model = model.to(device)

    # --- Dataset ---
    dataset = TokenizedGenomeDataset(
        config['manifest_path'],
        max_genomes=config.get('max_genomes'),
        min_genus_count=config.get('min_genus_count', 2),
    )
    loader = DataLoader(
        dataset, batch_size=config['batch_size'], shuffle=True,
        num_workers=0, drop_last=True,
    )

    # --- Optimizer: separate groups ---
    # Freeze encoder to isolate κ dynamics (optional)
    encoder_params = []
    for name, p in model.named_parameters():
        if name == 'hyper.c':
            continue
        if config.get('freeze_encoder', False):
            p.requires_grad = False
        else:
            encoder_params.append(p)

    param_groups = [
        {'params': [model.hyper.c], 'lr': config.get('c_lr', 1e-2), 'weight_decay': 0},
    ]
    if encoder_params and not config.get('freeze_encoder', False):
        param_groups.append({
            'params': encoder_params,
            'lr': config.get('lr', 1e-4),
            'weight_decay': 1e-4,
        })

    optimizer = torch.optim.Adam(param_groups)

    # --- Training loop ---
    history = []
    t0 = time.time()

    for epoch in range(config['epochs']):
        model.train()
        epoch_losses = {'total': 0, 'mlm': 0, 'dec': 0, 'hex': 0}
        n_batches = 0

        for tokens, labels in loader:
            tokens = tokens.to(device)
            labels = torch.tensor(labels, device=device) if not isinstance(labels, torch.Tensor) else labels.to(device)

            optimizer.zero_grad()

            # THE KEY FIX: pass tax_ids to activate InfoNCE
            loss, logs = model(tokens, tax_ids=labels)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            # Ensure c stays positive
            with torch.no_grad():
                model.hyper.c.data.clamp_(min=0.01)

            for k in epoch_losses:
                epoch_losses[k] += logs.get(k, logs.get('total', 0))
            n_batches += 1

        # Epoch summary
        avg = {k: v / max(n_batches, 1) for k, v in epoch_losses.items()}
        c_val = model.hyper.c.item()
        c_grad = model.hyper.c.grad.item() if model.hyper.c.grad is not None else 0.0

        record = {
            'epoch': epoch,
            'c': c_val,
            'c_grad': c_grad,
            **avg,
            'time': time.time() - t0,
        }
        history.append(record)

        if epoch % config.get('log_every', 1) == 0:
            print(f"Epoch {epoch:3d} | κ={c_val:.6f} | ∇κ={c_grad:+.6f} | "
                  f"loss={avg['total']:.4f} | mlm={avg['mlm']:.4f} | "
                  f"hex={avg['hex']:.4f} | t={time.time()-t0:.0f}s")

    return {
        'c_init': c_init,
        'seed': config['seed'],
        'c_final': model.hyper.c.item(),
        'loss_final': history[-1]['total'],
        'history': history,
        'n_genomes': len(dataset),
        'n_genera': dataset.n_genera,
        'freeze_encoder': config.get('freeze_encoder', False),
    }


# ─── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Measure κ via BiosphereCodec fine-tuning")
    parser.add_argument('--checkpoint', type=str, required=True)
    parser.add_argument('--codec-dir', type=str, required=True,
                        help='Directory containing BiosphereCodec.py')
    parser.add_argument('--manifest', type=str, required=True)
    parser.add_argument('--output', type=str, default='codec_kappa_results')
    parser.add_argument('--max-genomes', type=int, default=None)
    parser.add_argument('--epochs', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--d-model', type=int, default=512)
    parser.add_argument('--latent-dim', type=int, default=128)
    parser.add_argument('--vocab-size', type=int, default=4096)
    parser.add_argument('--n-layers', type=int, default=6)
    parser.add_argument('--max-len', type=int, default=131072)
    parser.add_argument('--lr', type=float, default=1e-4,
                        help='Encoder learning rate')
    parser.add_argument('--c-lr', type=float, default=1e-2,
                        help='Curvature learning rate')
    parser.add_argument('--c-inits', type=str, default='0.5,1.0,1.25,2.0,5.0')
    parser.add_argument('--seeds', type=str, default='42,59,76,93,110')
    parser.add_argument('--freeze-encoder', action='store_true',
                        help='Freeze encoder, only learn κ')
    parser.add_argument('--log-every', type=int, default=1)
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    c_inits = [float(x) for x in args.c_inits.split(',')]
    seeds = [int(x) for x in args.seeds.split(',')]

    base_config = {
        'checkpoint_path': args.checkpoint,
        'codec_dir': args.codec_dir,
        'manifest_path': args.manifest,
        'max_genomes': args.max_genomes,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'd_model': args.d_model,
        'latent_dim': args.latent_dim,
        'vocab_size': args.vocab_size,
        'n_layers': args.n_layers,
        'max_len': args.max_len,
        'lr': args.lr,
        'c_lr': args.c_lr,
        'freeze_encoder': args.freeze_encoder,
        'log_every': args.log_every,
    }

    print("=" * 60)
    print("κ MEASUREMENT — Pretrained BiosphereCodec + InfoNCE")
    print("=" * 60)
    print(f"Checkpoint: {args.checkpoint}")
    print(f"Manual Poincaré distance (no geoopt)")
    print(f"Encoder: {'FROZEN' if args.freeze_encoder else 'fine-tuning'}")
    print(f"c_inits: {c_inits}")
    print(f"seeds:   {seeds}")
    print()

    all_results = []

    for c_init in c_inits:
        for seed in seeds:
            config = {**base_config, 'c_init': c_init, 'seed': seed}
            try:
                result = run_finetune(config)
                all_results.append(result)
            except Exception as e:
                import traceback
                print(f"FAILED: c_init={c_init}, seed={seed}: {e}")
                traceback.print_exc()
                continue

            # Save incrementally
            with open(outdir / 'results.json', 'w') as f:
                json.dump({
                    'config': base_config,
                    'runs': [{k: v for k, v in r.items() if k != 'history'}
                             for r in all_results],
                }, f, indent=2)

            with open(outdir / f'c{c_init:.2f}_s{seed}_history.json', 'w') as f:
                json.dump(result['history'], f, indent=2)

    # ─── Summary ────────────────────────────────────────────────
    if all_results:
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        print(f"{'c_init':>8s} {'seed':>6s} {'c_final':>10s} {'loss':>8s} {'hex':>8s}")
        print(f"{'─'*8} {'─'*6} {'─'*10} {'─'*8} {'─'*8}")

        c_finals = []
        for r in all_results:
            hex_final = r['history'][-1].get('hex', 0) if r['history'] else 0
            print(f"{r['c_init']:>8.2f} {r['seed']:>6d} {r['c_final']:>10.6f} "
                  f"{r['loss_final']:>8.4f} {hex_final:>8.4f}")
            c_finals.append(r['c_final'])

        c_finals = np.array(c_finals)
        cv = c_finals.std() / c_finals.mean() * 100 if c_finals.mean() > 0 else 999
        print(f"\nκ = {c_finals.mean():.4f} ± {c_finals.std():.4f} (CV = {cv:.1f}%)")

        kappa_theory = (1.6 * math.log(2)) ** 2
        agreement = abs(c_finals.mean() - kappa_theory) / kappa_theory * 100
        print(f"Theory: κ = (1.6·ln2)² = {kappa_theory:.4f}")
        print(f"Agreement: {agreement:.1f}%")

        with open(outdir / 'summary.json', 'w') as f:
            json.dump({
                'c_final_mean': float(c_finals.mean()),
                'c_final_std': float(c_finals.std()),
                'c_final_cv': float(cv),
                'kappa_theory': float(kappa_theory),
                'agreement_pct': float(agreement),
                'runs': [{k: v for k, v in r.items() if k != 'history'}
                         for r in all_results],
            }, f, indent=2)


if __name__ == '__main__':
    main()
