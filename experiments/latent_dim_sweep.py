#!/usr/bin/env python3
"""
Latent Dimension Sweep Experiment
==================================

Separates two packing problems that CCS sometimes treated as one:

- Phylogenetic (HEX) structure is a tree metric. Its embeddability floor
  is n=2 (radial depth + angular divergence). HEX loss should plateau
  once the latent chart has that room; extra angular dimensions are
  reticulation / encoder slack, not more depth.
- Sequence reconstruction (MLM) is not a tree metric. It can demand a
  much larger latent dimension without touching the n=2 floor.

This sweep does **not** derive n=2 from the state equation, and it does
not treat a fitted kappa as a certified measurement (curvature is
degenerate with InfoNCE temperature; freeze it).

Experiment Design:
- Train BiosphereCodec with latent_dim ∈ {2, 4, 8, 16, 32, 64, 128}
- Two modes:
  1. HEX-only (phylogenetic structure only)
  2. HEX + MLM (full reconstruction)
- Measure:
  - HEX loss (phylogenetic fidelity) — expect plateau near dim 2–4
  - MLM loss (reconstruction fidelity) — may keep improving well above 2
  - If kappa is reported, treat it as a frozen design parameter, not a
    discovery

Prediction:
- HEX loss plateaus at dim ≈ 2-4 (embeddability, not a fitted n)
- MLM loss degrades sharply below dim ≈ 16-32 (different packing problem)
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

# Add parent to path for model imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from model.biosphere_codec import BiosphereCodec


class GenomicDataset(Dataset):
    """Simple dataset for genomic sequences."""

    def __init__(self, sequences: List[str], max_len: int = 1024):
        self.sequences = sequences
        self.max_len = max_len
        self.vocab = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5, '<pad>': 0}

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx][:self.max_len]
        tokens = [self.vocab.get(c, 5) for c in seq.upper()]
        # Pad
        if len(tokens) < self.max_len:
            tokens += [0] * (self.max_len - len(tokens))
        return torch.tensor(tokens, dtype=torch.long)


def train_one_epoch(
    model: BiosphereCodec,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    hex_only: bool = False,
) -> Dict[str, float]:
    """Train for one epoch, tracking losses separately."""
    model.train()

    total_loss = 0.0
    total_hex = 0.0
    total_mlm = 0.0
    n_batches = 0

    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()

        # Forward pass
        outputs = model(batch)

        if isinstance(outputs, tuple):
            loss, losses_dict = outputs
        else:
            loss = outputs
            losses_dict = {}

        # Track individual losses
        hex_loss = losses_dict.get('hex', 0.0)
        mlm_loss = losses_dict.get('mlm', 0.0)

        if hex_only:
            # Only backprop through HEX loss
            if 'hex' in losses_dict:
                loss = losses_dict['hex']
            else:
                loss = loss  # fallback

        loss.backward()
        optimizer.step()

        total_loss += loss.item()
        total_hex += hex_loss if isinstance(hex_loss, float) else hex_loss.item() if hasattr(hex_loss, 'item') else 0
        total_mlm += mlm_loss if isinstance(mlm_loss, float) else mlm_loss.item() if hasattr(mlm_loss, 'item') else 0
        n_batches += 1

    return {
        'loss': total_loss / max(n_batches, 1),
        'hex': total_hex / max(n_batches, 1),
        'mlm': total_mlm / max(n_batches, 1),
    }


def extract_kappa(model: BiosphereCodec) -> float:
    """Extract learned curvature from model."""
    # Try different attribute names
    for attr_path in ['hyper.c', 'hyper.manifold.c', 'hyperbolic.c', 'poincare.c']:
        try:
            obj = model
            for part in attr_path.split('.'):
                obj = getattr(obj, part)
            if hasattr(obj, 'item'):
                return obj.item()
            return float(obj)
        except AttributeError:
            continue
    return float('nan')


def run_experiment(
    latent_dim: int,
    sequences: List[str],
    n_epochs: int = 50,
    batch_size: int = 4,
    hex_only: bool = False,
    device: str = 'cpu',
) -> Dict:
    """Run single experiment with given latent dimension."""

    print(f"\n{'='*60}")
    print(f"Running: latent_dim={latent_dim}, hex_only={hex_only}")
    print(f"{'='*60}")

    # Create dataset
    dataset = GenomicDataset(sequences, max_len=1024)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    # Create model
    model = BiosphereCodec(
        vocab=6,
        d_model=64,  # Smaller for faster experimentation
        n_layers=2,
        latent_dim=latent_dim,
    )
    model = model.to(device)

    # Optimizer
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)

    # Training loop
    history = {
        'epoch': [],
        'loss': [],
        'hex': [],
        'mlm': [],
        'kappa': [],
    }

    for epoch in range(n_epochs):
        metrics = train_one_epoch(model, loader, optimizer, device, hex_only)
        kappa = extract_kappa(model)

        history['epoch'].append(epoch)
        history['loss'].append(metrics['loss'])
        history['hex'].append(metrics['hex'])
        history['mlm'].append(metrics['mlm'])
        history['kappa'].append(kappa)

        if epoch % 10 == 0 or epoch == n_epochs - 1:
            print(f"  Epoch {epoch:3d}: loss={metrics['loss']:.4f}, "
                  f"hex={metrics['hex']:.4f}, mlm={metrics['mlm']:.4f}, κ={kappa:.4f}")

    # Final metrics
    return {
        'latent_dim': latent_dim,
        'hex_only': hex_only,
        'final_loss': history['loss'][-1],
        'final_hex': history['hex'][-1],
        'final_mlm': history['mlm'][-1],
        'final_kappa': history['kappa'][-1],
        'kappa_trajectory': history['kappa'],
        'history': history,
    }


def main():
    parser = argparse.ArgumentParser(description='Latent Dimension Sweep Experiment')
    parser.add_argument('--fasta', type=str, help='Path to FASTA file with sequences')
    parser.add_argument('--dims', type=str, default='2,4,8,16,32,64,128',
                        help='Comma-separated latent dimensions to test')
    parser.add_argument('--epochs', type=int, default=50, help='Epochs per experiment')
    parser.add_argument('--batch-size', type=int, default=4, help='Batch size')
    parser.add_argument('--hex-only', action='store_true', help='Train with HEX loss only')
    parser.add_argument('--device', type=str, default='cpu', help='Device (cpu/cuda)')
    parser.add_argument('--output', type=str, default='latent_dim_results.json',
                        help='Output file for results')
    args = parser.parse_args()

    # Load sequences
    if args.fasta and Path(args.fasta).exists():
        from Bio import SeqIO
        sequences = [str(rec.seq) for rec in SeqIO.parse(args.fasta, 'fasta')]
        print(f"Loaded {len(sequences)} sequences from {args.fasta}")
    else:
        # Generate synthetic sequences for testing
        print("No FASTA provided, using synthetic sequences...")
        np.random.seed(42)
        sequences = [
            ''.join(np.random.choice(list('ACGT'), size=2000))
            for _ in range(100)
        ]

    # Parse dimensions
    dims = [int(d) for d in args.dims.split(',')]

    # Run experiments
    results = []
    for dim in dims:
        result = run_experiment(
            latent_dim=dim,
            sequences=sequences,
            n_epochs=args.epochs,
            batch_size=args.batch_size,
            hex_only=args.hex_only,
            device=args.device,
        )
        results.append(result)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'Dim':>6} {'Loss':>10} {'HEX':>10} {'MLM':>10} {'κ':>10}")
    print("-"*50)
    for r in results:
        print(f"{r['latent_dim']:>6} {r['final_loss']:>10.4f} {r['final_hex']:>10.4f} "
              f"{r['final_mlm']:>10.4f} {r['final_kappa']:>10.4f}")

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert numpy arrays to lists for JSON serialization
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        return obj

    results_json = json.loads(json.dumps(results, default=convert))

    with open(output_path, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'args': vars(args),
            'results': results_json,
        }, f, indent=2)

    print(f"\nResults saved to: {output_path}")


if __name__ == '__main__':
    main()
