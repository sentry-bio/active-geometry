#!/usr/bin/env python3
"""
Coordinate Convergence — Real Procrustes Analysis on Actual Checkpoints
=======================================================================

This script loads ALL 5 trained checkpoints, encodes the SAME set of genomes
through each, and computes pairwise Procrustes correlations to measure
coordinate convergence.

The 5-seed training runs (June-July 2025) used FIXED κ=1.0 (InfoNCE/DIST
losses omitted), so this measures COORDINATE convergence only, not curvature.

Seeds: original(0), 42, 137, 2024, 888
"""

import os
import sys
import json
import pickle
import argparse
import logging
import time
from pathlib import Path
from typing import Dict, List, Tuple
from itertools import combinations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
import numpy as np
import zstandard as zstd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
log = logging.getLogger("Procrustes")


# =========================================================================
# 1. MODEL (must match training architecture exactly)
# =========================================================================

class HyenaOperator(nn.Module):
    def __init__(self, d_model: int, mode: str = "encoder"):
        super().__init__()
        # Conv1d with bias=False (checkpoint has no depthwise.bias)
        self.depthwise = nn.Conv1d(d_model, d_model, 7, padding=3, groups=d_model, bias=False)
        # Gate is Conv1d kernel=1, not Linear (checkpoint shape: [256, 256, 1])
        self.gate = nn.Conv1d(d_model, d_model, 1)
        self.mode = mode

    def forward(self, x: Tensor) -> Tensor:
        # x: [B, L, D]
        xt = x.transpose(1, 2)  # [B, D, L]
        gate = torch.sigmoid(self.gate(xt)).transpose(1, 2)  # [B, L, D]
        conv = self.depthwise(xt).transpose(1, 2)  # [B, L, D]
        return gate * conv


class AttentionPool(nn.Module):
    """Matches checkpoint key: encoder.pool.attn_gate"""
    def __init__(self, d_model: int):
        super().__init__()
        self.attn_gate = nn.Linear(d_model, 1)

    def forward(self, h: Tensor) -> Tensor:
        mean_pool = h.mean(dim=1)
        max_pool = h.max(dim=1).values
        attn_w = torch.softmax(self.attn_gate(h).squeeze(-1), dim=1)
        attn_pool = (h * attn_w.unsqueeze(-1)).sum(dim=1)
        return torch.cat([mean_pool, max_pool, attn_pool], dim=-1)


class BiosphereEncoder(nn.Module):
    def __init__(self, vocab: int = 5444, d_model: int = 256,
                 n_layers: int = 4, max_len: int = 8192):
        super().__init__()
        self.embed = nn.Embedding(vocab, d_model)
        self.pos = nn.Parameter(torch.randn(max_len, d_model) * 0.02)

        self.layers = nn.ModuleList()
        for _ in range(n_layers):
            self.layers.append(nn.ModuleDict({
                'norm1': nn.LayerNorm(d_model),
                'hyena': HyenaOperator(d_model),
                'norm2': nn.LayerNorm(d_model),
                'ff': nn.Sequential(
                    nn.Linear(d_model, d_model * 4),
                    nn.GELU(), nn.Dropout(0.0),  # No dropout at inference
                    nn.Linear(d_model * 4, d_model),
                ),
            }))

        self.norm = nn.LayerNorm(d_model)
        self.pool = AttentionPool(d_model)

    def forward(self, ids: Tensor) -> Tuple[Tensor, Tensor]:
        B, L = ids.shape
        x = self.embed(ids) + self.pos[:L]
        for layer in self.layers:
            x = x + layer['hyena'](layer['norm1'](x))
            x = x + layer['ff'](layer['norm2'](x))
        h = self.norm(x)
        pooled = self.pool(h)
        return h, pooled


class PoincareMapping(nn.Module):
    def __init__(self, in_dim: int = 768, latent_dim: int = 256, init_c: float = 1.0):
        super().__init__()
        self.lin = nn.Linear(in_dim, latent_dim)
        self.c = nn.Parameter(torch.tensor(init_c))

    def forward(self, x: Tensor) -> Tensor:
        z = torch.tanh(self.lin(x))
        c = torch.clamp(self.c, min=1e-4)
        r_max = 0.9 / torch.sqrt(c)
        z = z * r_max
        norm = torch.norm(z, dim=-1, keepdim=True).clamp(min=1e-8)
        max_norm = r_max - 1e-5
        factor = torch.where(norm > max_norm, max_norm / norm, torch.ones_like(norm))
        return z * factor


class BiosphereDecoder(nn.Module):
    """Lightweight decoder matching checkpoint keys."""
    def __init__(self, shared_embed: nn.Embedding):
        super().__init__()
        d_model = shared_embed.embedding_dim
        self.hyena = HyenaOperator(d_model, mode="causal")
        self.norm = nn.LayerNorm(d_model)
        self.proj = nn.Linear(d_model, shared_embed.num_embeddings, bias=False)
        self.proj.weight = shared_embed.weight

    def forward(self, h: Tensor) -> Tensor:
        return self.proj(self.norm(self.hyena(h)))


class LossFnModule(nn.Module):
    """Placeholder for loss_fn.manifold in checkpoint."""
    def __init__(self):
        super().__init__()
        self.manifold = PoincareMapping()


class BiosphereCodecInference(nn.Module):
    """Inference-only model matching the 5-seed training checkpoints."""
    def __init__(self):
        super().__init__()
        self.encoder = BiosphereEncoder()
        self.hyper = PoincareMapping()
        self.decoder = BiosphereDecoder(self.encoder.embed)
        self.loss_fn = LossFnModule()

    def encode(self, ids: Tensor) -> Tensor:
        """Encode token IDs → Poincaré ball coordinates."""
        _, pooled = self.encoder(ids)
        z = self.hyper(pooled)
        return z


# =========================================================================
# 2. DATA LOADING
# =========================================================================

def load_tokenizer(path: str):
    with open(path, 'rb') as f:
        tok = pickle.load(f)
    return tok['vocabulary']


def tokenize_sequence(seq: str, vocab: dict, max_len: int = 8192) -> torch.Tensor:
    pad_id = vocab.get('<PAD>', 0)
    unk_id = vocab.get('<UNK>', 1)
    tokens = []
    for i in range(len(seq) - 2):
        kmer = seq[i:i+3]
        tokens.append(vocab.get(kmer, unk_id))
    if len(tokens) > max_len:
        tokens = tokens[:max_len]
    elif len(tokens) < max_len:
        tokens = tokens + [pad_id] * (max_len - len(tokens))
    return torch.tensor(tokens, dtype=torch.long)


def load_sequences(data_root: str, max_samples: int = 200):
    """Load a fixed subset of sequences for comparison across models."""
    sequences = []
    genome_ids = []
    data_path = Path(data_root)
    for zst_file in sorted(data_path.glob("**/*.zst")):
        if len(sequences) >= max_samples:
            break
        with open(zst_file, 'rb') as f:
            dctx = zstd.ZstdDecompressor()
            raw = dctx.decompress(f.read())
            data = json.loads(raw.decode('utf-8'))
        for item in data:
            if len(sequences) >= max_samples:
                break
            if 'sequence' not in item or len(item['sequence']) < 100:
                continue
            sequences.append(item['sequence'][:8192])
            genome_ids.append(item.get('genome_id', f'seq_{len(sequences)}'))
    return sequences, genome_ids


# =========================================================================
# 3. PROCRUSTES ANALYSIS
# =========================================================================

def procrustes(X: np.ndarray, Y: np.ndarray) -> Tuple[float, np.ndarray]:
    """Ordinary Procrustes analysis.

    Returns:
        correlation: Procrustes correlation (1 - normalized disparity)
        Y_aligned: Y rotated/reflected to best match X
    """
    # Center
    X_c = X - X.mean(axis=0)
    Y_c = Y - Y.mean(axis=0)

    # Scale
    X_s = np.sqrt((X_c ** 2).sum())
    Y_s = np.sqrt((Y_c ** 2).sum())
    X_c = X_c / X_s
    Y_c = Y_c / Y_s

    # Optimal rotation via SVD
    M = X_c.T @ Y_c
    U, S, Vt = np.linalg.svd(M)
    R = Vt.T @ U.T

    # Correlation from singular values
    correlation = S.sum()

    Y_aligned = Y_c @ R
    return correlation, Y_aligned


# =========================================================================
# 4. MAIN
# =========================================================================

def main():
    parser = argparse.ArgumentParser(description="Real Procrustes analysis on 5-seed checkpoints")
    parser.add_argument("--data_root", type=str,
                        default="/zfs_raid/SentryBio/5k_test_genomes/processed_biosphere_data_supervised")
    parser.add_argument("--tokenizer_path", type=str,
                        default="/zfs_raid/SentryBio/5k_test_genomes/tokenizer_output_bulletproof/biosphere_tokenizer_bulletproof.pkl")
    parser.add_argument("--n_samples", type=int, default=200,
                        help="Number of genomes to embed for comparison")
    parser.add_argument("--output_dir", type=str, default="./procrustes_results")
    parser.add_argument("--batch_size", type=int, default=4)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Device: {device}")

    # Checkpoint paths (THE REAL ONES)
    checkpoints = {
        'original': '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_final/checkpoint_step_7000.pt',
        'seed_42': '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_42/checkpoint_step_7000.pt',
        'seed_137': '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_137/checkpoint_step_7000.pt',
        'seed_2024': '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_2024/checkpoint_step_7000.pt',
        'seed_888': '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_888/checkpoint_step_7000.pt',
    }

    # Verify all checkpoints exist
    for name, path in checkpoints.items():
        if not Path(path).exists():
            log.error(f"MISSING checkpoint: {name} at {path}")
            sys.exit(1)
    log.info(f"All {len(checkpoints)} checkpoints verified")

    # Load tokenizer and sequences
    log.info("Loading tokenizer and sequences...")
    vocab = load_tokenizer(args.tokenizer_path)
    sequences, genome_ids = load_sequences(args.data_root, args.n_samples)
    log.info(f"Loaded {len(sequences)} sequences")

    # Tokenize all sequences once
    log.info("Tokenizing...")
    all_tokens = torch.stack([tokenize_sequence(seq, vocab) for seq in sequences])
    log.info(f"Token tensor: {all_tokens.shape}")

    # Encode with each model
    embeddings = {}
    for name, ckpt_path in checkpoints.items():
        log.info(f"Encoding with {name}...")
        model = BiosphereCodecInference().to(device)
        ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(ckpt['model_state_dict'], strict=True)
        model.eval()

        all_z = []
        with torch.no_grad():
            for i in range(0, len(all_tokens), args.batch_size):
                batch = all_tokens[i:i+args.batch_size].to(device)
                z = model.encode(batch)
                all_z.append(z.cpu().numpy())

        embeddings[name] = np.concatenate(all_z, axis=0)
        kappa = ckpt['model_state_dict']['hyper.c'].item()
        log.info(f"  {name}: {embeddings[name].shape}, κ={kappa:.6f}")

        del model
        torch.cuda.empty_cache()

    # Pairwise Procrustes
    log.info("\n" + "="*60)
    log.info("PAIRWISE PROCRUSTES ANALYSIS")
    log.info("="*60)

    names = list(checkpoints.keys())
    correlations = {}
    for i, j in combinations(range(len(names)), 2):
        n1, n2 = names[i], names[j]
        corr, _ = procrustes(embeddings[n1], embeddings[n2])
        correlations[f"{n1}_vs_{n2}"] = corr
        log.info(f"  {n1:>12s} vs {n2:<12s}: r = {corr:.6f}")

    all_corrs = list(correlations.values())
    mean_corr = np.mean(all_corrs)
    std_corr = np.std(all_corrs)
    min_corr = np.min(all_corrs)
    max_corr = np.max(all_corrs)

    log.info(f"\nSUMMARY:")
    log.info(f"  Mean Procrustes correlation: {mean_corr:.6f} ± {std_corr:.6f}")
    log.info(f"  Range: [{min_corr:.6f}, {max_corr:.6f}]")
    log.info(f"  All > 0.99: {all(c > 0.99 for c in all_corrs)}")
    log.info(f"  All > 0.95: {all(c > 0.95 for c in all_corrs)}")
    log.info(f"  All > 0.90: {all(c > 0.90 for c in all_corrs)}")

    # Also compute cosine similarity of mean embedding directions
    log.info("\nMEAN EMBEDDING COSINE SIMILARITY:")
    for i, j in combinations(range(len(names)), 2):
        n1, n2 = names[i], names[j]
        m1 = embeddings[n1].mean(axis=0)
        m2 = embeddings[n2].mean(axis=0)
        cos = np.dot(m1, m2) / (np.linalg.norm(m1) * np.linalg.norm(m2) + 1e-10)
        log.info(f"  {n1:>12s} vs {n2:<12s}: cos = {cos:.6f}")

    # Embedding statistics
    log.info("\nEMBEDDING STATISTICS:")
    for name in names:
        z = embeddings[name]
        norms = np.linalg.norm(z, axis=1)
        log.info(f"  {name:>12s}: mean_norm={norms.mean():.4f} ± {norms.std():.4f}, "
                 f"range=[{norms.min():.4f}, {norms.max():.4f}]")

    # Save results
    # Convert numpy floats to Python floats for JSON
    correlations_clean = {k: float(v) for k, v in correlations.items()}

    results = {
        'experiment': 'Coordinate Convergence — Real Procrustes Analysis',
        'date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'n_sequences': len(sequences),
        'n_models': len(checkpoints),
        'seeds': {name: path for name, path in checkpoints.items()},
        'kappa_values': {name: 1.0 for name in names},  # All fixed at 1.0
        'pairwise_procrustes': correlations_clean,
        'summary': {
            'mean_correlation': float(mean_corr),
            'std_correlation': float(std_corr),
            'min_correlation': float(min_corr),
            'max_correlation': float(max_corr),
            'all_above_099': bool(all(c > 0.99 for c in all_corrs)),
            'all_above_095': bool(all(c > 0.95 for c in all_corrs)),
            'all_above_090': bool(all(c > 0.90 for c in all_corrs)),
        },
        'genome_ids': genome_ids[:20],  # First 20 for reference
        'note': 'All 5 checkpoints have kappa=1.0 (fixed). '
                'This measures COORDINATE convergence, not curvature convergence.',
    }

    results_path = output_dir / "procrustes_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    log.info(f"\nResults saved to {results_path}")

    # Save embedding matrices for future analysis
    np.savez(output_dir / "embeddings.npz",
             **{name: emb for name, emb in embeddings.items()})
    log.info(f"Embeddings saved to {output_dir / 'embeddings.npz'}")


if __name__ == "__main__":
    main()
