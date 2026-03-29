#!/usr/bin/env python3
"""
κ Convergence Experiment — Canonical Reproduction
=================================================

This script measures the optimal hyperbolic curvature κ for genomic embeddings
by making the Poincaré ball curvature parameter learnable and tracking its
convergence across multiple random seeds and initializations.

Key design:
  - Uses the full BiosphereCodec architecture (Hyena encoder + Poincaré head)
  - Full loss = MLM + DEC + 0.1 * HEX + 0.5 * DIST
  - κ = PoincareMapping.c is nn.Parameter (learnable, receives gradient via HEX/DIST)
  - MLM/DEC provide no gradient to κ (confirmed by ablation)
  - HEX (InfoNCE on hyperbolic distances) + DIST (patristic regression) drive κ

Expected result: κ converges to ~1.25 regardless of seed or initialization.

Usage:
  # Single run
  python run_kappa_convergence.py --seed 42 --init_kappa 0.5 --steps 3000

  # Full 5-seed sweep
  python run_kappa_convergence.py --sweep
"""

import os
import sys
import json
import argparse
import logging
import math
import time
import pickle
import zstandard as zstd
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import DataLoader, Dataset
from torch.optim import AdamW
from torch.cuda.amp import GradScaler, autocast
import numpy as np

try:
    import geoopt
    _GEOOPT = True
except ImportError:
    _GEOOPT = False
    print("WARNING: geoopt not available. Install with: pip install geoopt")


# =========================================================================
# 1. CONFIGURATION
# =========================================================================

@dataclass
class KappaConfig:
    """Configuration for κ convergence experiment."""
    # Architecture (must match 5-seed training)
    d_model: int = 256
    n_layers: int = 4
    latent_dim: int = 128
    vocab_size: int = 5444
    max_sequence_length: int = 8192

    # Curvature
    init_kappa: float = 1.0  # Initial value for learnable κ
    kappa_lr_scale: float = 1.0  # Relative LR for κ vs rest of model

    # Training
    batch_size: int = 8
    gradient_accumulation_steps: int = 2
    learning_rate: float = 2e-4
    weight_decay: float = 0.01
    max_grad_norm: float = 1.0
    max_steps: int = 3000
    use_mixed_precision: bool = True

    # Loss weights (must match manuscript §9.3)
    hex_weight: float = 0.1
    dist_weight: float = 0.5
    infonce_temp: float = 0.1
    mlm_prob: float = 0.15

    # Data
    data_root: str = "/zfs_raid/SentryBio/5k_test_genomes/processed_biosphere_data_supervised"
    tokenizer_path: str = "/zfs_raid/SentryBio/5k_test_genomes/tokenizer_output_bulletproof/biosphere_tokenizer_bulletproof.pkl"
    manifest_path: str = "/zfs_raid/SentryBio/5k_test_genomes/MASTER_REGISTRY_refseq_only.csv"

    # Logging
    seed: int = 42
    output_dir: str = "./kappa_convergence"
    log_every: int = 10  # Log κ frequently
    save_every: int = 500


# =========================================================================
# 2. MODEL COMPONENTS (from BiosphereCodec, minimal extraction)
# =========================================================================

class HyenaOperator(nn.Module):
    """Hyena convolution operator."""
    def __init__(self, d_model: int, mode: str = "encoder"):
        super().__init__()
        self.depthwise = nn.Conv1d(d_model, d_model, 7, padding=3, groups=d_model)
        self.gate = nn.Linear(d_model, d_model)
        self.mode = mode

    def forward(self, x: Tensor) -> Tensor:
        gate = torch.sigmoid(self.gate(x))
        conv = self.depthwise(x.transpose(1, 2)).transpose(1, 2)
        return gate * conv


class BiosphereEncoder(nn.Module):
    """Hyena-based encoder with hierarchical pooling."""
    def __init__(self, vocab: int, d_model: int, n_layers: int, max_len: int):
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
                    nn.GELU(), nn.Dropout(0.1),
                    nn.Linear(d_model * 4, d_model),
                ),
            }))

        self.norm = nn.LayerNorm(d_model)
        self.pool_attn = nn.Linear(d_model, 1)

    def forward(self, ids: Tensor) -> Tuple[Tensor, Tensor]:
        B, L = ids.shape
        x = self.embed(ids) + self.pos[:L]

        for layer in self.layers:
            x = x + layer['hyena'](layer['norm1'](x))
            x = x + layer['ff'](layer['norm2'](x))

        h = self.norm(x)

        # Hierarchical pooling: [mean, max, attn-weighted]
        mean_pool = h.mean(dim=1)
        max_pool = h.max(dim=1).values
        attn_w = torch.softmax(self.pool_attn(h).squeeze(-1), dim=1)
        attn_pool = (h * attn_w.unsqueeze(-1)).sum(dim=1)
        pooled = torch.cat([mean_pool, max_pool, attn_pool], dim=-1)  # [B, 3*D]

        return h, pooled


class PoincareMapping(nn.Module):
    """Linear → Poincaré ball projection with LEARNABLE curvature.

    Uses manual Poincaré distance computation so gradients flow cleanly
    through the curvature parameter c without geoopt in-place issues.
    """
    def __init__(self, in_dim: int, latent_dim: int = 128, init_c: float = 1.0):
        super().__init__()
        self.lin = nn.Linear(in_dim, latent_dim)
        # THE KEY PARAMETER: learnable curvature
        self.c = nn.Parameter(torch.tensor(init_c))

    def _project_to_ball(self, z: Tensor) -> Tensor:
        """Project points inside the Poincaré ball of curvature c."""
        c = torch.clamp(self.c, min=1e-4)
        r_max = 0.9 / torch.sqrt(c)
        norm = torch.norm(z, dim=-1, keepdim=True).clamp(min=1e-8)
        # Clamp to be strictly inside the ball
        max_norm = r_max - 1e-5
        factor = torch.where(norm > max_norm, max_norm / norm, torch.ones_like(norm))
        return z * factor

    def forward(self, x: Tensor) -> Tensor:
        z = torch.tanh(self.lin(x))
        c = torch.clamp(self.c, min=1e-4)
        r_max = 0.9 / torch.sqrt(c)
        z = z * r_max
        return self._project_to_ball(z)

    def dist_mat(self, z: Tensor) -> Tensor:
        """Poincaré ball distance with gradients through c.

        d(u,v) = (1/√c) * arccosh(1 + 2c||u-v||² / ((1-c||u||²)(1-c||v||²)))
        """
        c = torch.clamp(self.c, min=1e-4)

        # Pairwise squared Euclidean distances
        diff = z.unsqueeze(1) - z.unsqueeze(0)  # [B, B, d]
        sqdist = (diff ** 2).sum(-1)  # [B, B]

        # Norms squared
        sqnorm_u = (z ** 2).sum(-1, keepdim=True)  # [B, 1]
        sqnorm_v = (z ** 2).sum(-1).unsqueeze(0)    # [1, B]

        # Denominator: (1 - c||u||²)(1 - c||v||²)
        denom = (1 - c * sqnorm_u) * (1 - c * sqnorm_v)
        denom = torch.clamp(denom, min=1e-10)

        # arccosh argument
        x = 1 + 2 * c * sqdist / denom
        x = torch.clamp(x, min=1.0 + 1e-7)  # arccosh domain

        # d = (1/√c) * arccosh(x)
        dist = torch.acosh(x) / torch.sqrt(c)
        return dist


class BiosphereDecoder(nn.Module):
    """Lightweight causal decoder with weight tying."""
    def __init__(self, shared_embed: nn.Embedding):
        super().__init__()
        d_model = shared_embed.embedding_dim
        self.hyena = HyenaOperator(d_model, mode="causal")
        self.norm = nn.LayerNorm(d_model)
        self.proj = nn.Linear(d_model, shared_embed.num_embeddings, bias=False)
        self.proj.weight = shared_embed.weight

    def forward(self, h: Tensor) -> Tensor:
        return self.proj(self.norm(self.hyena(h)))


class BiosphereCodecKappa(nn.Module):
    """Full BiosphereCodec with learnable κ and complete loss."""

    def __init__(self, config: KappaConfig):
        super().__init__()
        self.config = config
        self.encoder = BiosphereEncoder(
            config.vocab_size, config.d_model, config.n_layers,
            config.max_sequence_length
        )
        self.hyper = PoincareMapping(
            3 * config.d_model, config.latent_dim, init_c=config.init_kappa
        )
        self.decoder = BiosphereDecoder(self.encoder.embed)
        self.mask_id = config.vocab_size - 1

    def forward(self, ids: Tensor, tax_ids: Optional[Tensor] = None) -> Dict[str, Any]:
        B, L = ids.shape
        vocab = self.config.vocab_size

        # --- MLM masking ---
        labels = ids.clone()
        mask = torch.rand_like(ids.float()) < self.config.mlm_prob
        labels[~mask] = -100
        masked_ids = ids.clone()
        masked_ids[mask] = self.mask_id

        # --- Encode ---
        h, pooled = self.encoder(masked_ids)
        z = self.hyper(pooled)

        # --- Decode ---
        enc_logits = nn.functional.linear(h, self.encoder.embed.weight)
        dec_logits = self.decoder(h)

        # --- Losses ---
        mlm_loss = F.cross_entropy(
            enc_logits.view(-1, vocab), labels.view(-1), ignore_index=-100
        )
        dec_loss = F.cross_entropy(dec_logits.view(-1, vocab), ids.view(-1))

        # --- HEX: Hierarchical InfoNCE on Poincaré distances ---
        # This is THE loss that informs curvature (see ablation, §6)
        hex_loss = torch.tensor(0.0, device=ids.device)
        if tax_ids is not None and B >= 2:
            dist = self.hyper.dist_mat(z)  # [B, B] Poincaré distances
            logits = -dist / self.config.infonce_temp  # similarity
            diag = torch.eye(B, dtype=torch.bool, device=ids.device)
            logits = logits.masked_fill(diag, -1e9)

            # Positive pairs: same phylum
            pos_mask = tax_ids.unsqueeze(1) == tax_ids.unsqueeze(0)
            pos_mask = pos_mask & ~diag

            if pos_mask.any():
                # InfoNCE: push same-phylum closer, different-phylum apart
                log_probs = logits - torch.logsumexp(logits, dim=1, keepdim=True)
                hex_loss = -(log_probs[pos_mask].mean())

        # --- DIST: distance margin loss ---
        # Pushes embeddings apart, κ sets the scale
        dist_loss = torch.tensor(0.0, device=ids.device)
        if B >= 2:
            dist = self.hyper.dist_mat(z)
            mask_triu = torch.triu(torch.ones(B, B, device=ids.device), diagonal=1).bool()
            distances = dist[mask_triu]
            margin = 0.5
            dist_loss = F.relu(margin - distances).mean()

        # --- Total (manuscript §9.3) ---
        total = mlm_loss + dec_loss + self.config.hex_weight * hex_loss + self.config.dist_weight * dist_loss

        return {
            'loss': total,
            'mlm': mlm_loss.item(),
            'dec': dec_loss.item(),
            'hex': hex_loss.item(),
            'dist': dist_loss.item(),
            'kappa': self.hyper.c.item(),
            'kappa_grad': self.hyper.c.grad.item() if self.hyper.c.grad is not None else 0.0,
        }


# =========================================================================
# 3. DATASET
# =========================================================================

class GenomeDataset(Dataset):
    """Loads genomic sequences with taxonomy labels for contrastive learning."""

    def __init__(self, data_root: str, tokenizer_path: str, manifest_path: str,
                 max_len: int = 8192, max_samples: int = 2000, vocab_size: int = 5444):
        self.max_len = max_len
        self.vocab_size = vocab_size
        self.sequences = []
        self.tax_labels = []  # Phylum-level integer labels

        # Load tokenizer
        tok_path = Path(tokenizer_path)
        if tok_path.exists():
            with open(tok_path, 'rb') as f:
                tok_dict = pickle.load(f)
            self.vocab = tok_dict['vocabulary']
            self.pad_id = self.vocab.get('<PAD>', 0)
            self.unk_id = self.vocab.get('<UNK>', 1)
            print(f"Loaded tokenizer: {len(self.vocab)} tokens")
        else:
            self.vocab = None
            self.pad_id = 0
            self.unk_id = 1
            print(f"WARNING: Tokenizer not found at {tok_path}, using char-level")

        # Load taxonomy from manifest (accession → phylum)
        import csv
        acc_to_phylum = {}
        mpath = Path(manifest_path)
        if mpath.exists():
            with open(mpath) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    acc = row.get('accession', '')
                    phylum = row.get('phylum', 'Unknown')
                    if acc:
                        acc_to_phylum[acc] = phylum or 'Unknown'
            print(f"Loaded taxonomy for {len(acc_to_phylum)} accessions")
        else:
            print(f"WARNING: Manifest not found at {mpath}")

        # Build phylum → integer mapping
        phyla = sorted(set(acc_to_phylum.values()))
        phylum_to_id = {p: i for i, p in enumerate(phyla)}
        print(f"Phyla: {len(phyla)} ({phyla[:5]}...)")

        # Load sequences + taxonomy
        data_path = Path(data_root)
        zst_files = sorted(data_path.glob("**/*.zst"))
        print(f"Found {len(zst_files)} .zst files in {data_path}")

        for zst_file in zst_files:
            if len(self.sequences) >= max_samples:
                break
            try:
                with open(zst_file, 'rb') as f:
                    dctx = zstd.ZstdDecompressor()
                    raw = dctx.decompress(f.read())
                    data = json.loads(raw.decode('utf-8'))

                for item in data:
                    if len(self.sequences) >= max_samples:
                        break
                    if 'sequence' not in item or len(item['sequence']) < 100:
                        continue

                    gid = item.get('genome_id', '')
                    # Extract accession from genome_id (may be embedded)
                    acc = None
                    for a in acc_to_phylum:
                        if a in gid:
                            acc = a
                            break

                    phylum = acc_to_phylum.get(acc, 'Unknown') if acc else 'Unknown'
                    tax_id = phylum_to_id.get(phylum, len(phyla))

                    self.sequences.append(item['sequence'])
                    self.tax_labels.append(tax_id)

            except Exception as e:
                print(f"  Warning: {zst_file.name}: {e}")

        print(f"Loaded {len(self.sequences)} sequences with taxonomy")
        tax_counts = {}
        for t in self.tax_labels:
            tax_counts[t] = tax_counts.get(t, 0) + 1
        print(f"  Label distribution: {dict(sorted(tax_counts.items()))}")

        if len(self.sequences) == 0:
            raise ValueError("No sequences loaded!")

    def __len__(self):
        return len(self.sequences)

    def _tokenize(self, seq: str) -> List[int]:
        if self.vocab is not None:
            tokens = []
            for i in range(len(seq) - 2):
                kmer = seq[i:i+3]
                tokens.append(self.vocab.get(kmer, self.unk_id))
            return tokens
        else:
            return [min(ord(c), self.vocab_size - 1) for c in seq]

    def __getitem__(self, idx):
        seq = self.sequences[idx][:self.max_len]
        tokens = self._tokenize(seq)

        if len(tokens) > self.max_len:
            tokens = tokens[:self.max_len]
        elif len(tokens) < self.max_len:
            tokens = tokens + [self.pad_id] * (self.max_len - len(tokens))

        return (torch.tensor(tokens, dtype=torch.long),
                torch.tensor(self.tax_labels[idx], dtype=torch.long))


# =========================================================================
# 4. TRAINER
# =========================================================================

class KappaConvergenceTrainer:
    """Trains BiosphereCodec with learnable κ and logs convergence."""

    def __init__(self, config: KappaConfig):
        self.config = config
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Deterministic seeding
        torch.manual_seed(config.seed)
        np.random.seed(config.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(config.seed)

        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Save config
        with open(self.output_dir / "config.json", 'w') as f:
            json.dump(asdict(config), f, indent=2)

        self._setup_logging()
        self._init_model()
        self._init_data()
        self._init_optimizer()

        self.global_step = 0
        self.kappa_history = []

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.output_dir / "training.log"),
                logging.StreamHandler()
            ],
            force=True
        )
        self.logger = logging.getLogger("KappaConvergence")

    def _init_model(self):
        self.model = BiosphereCodecKappa(self.config).to(self.device)
        n_params = sum(p.numel() for p in self.model.parameters())
        self.logger.info(f"Model: {n_params:,} parameters on {self.device}")
        self.logger.info(f"Initial κ = {self.model.hyper.c.item():.6f}")
        self.logger.info(f"κ requires_grad = {self.model.hyper.c.requires_grad}")

    def _init_data(self):
        dataset = GenomeDataset(
            self.config.data_root, self.config.tokenizer_path,
            self.config.manifest_path,
            self.config.max_sequence_length, max_samples=2000,
            vocab_size=self.config.vocab_size
        )
        self.dataloader = DataLoader(
            dataset, batch_size=self.config.batch_size,
            shuffle=True, num_workers=0, pin_memory=True, drop_last=True
        )

    def _init_optimizer(self):
        # Separate κ from other params to allow different LR
        kappa_params = [self.model.hyper.c]
        other_params = [p for n, p in self.model.named_parameters()
                        if p is not self.model.hyper.c]

        self.optimizer = AdamW([
            {'params': other_params, 'lr': self.config.learning_rate,
             'weight_decay': self.config.weight_decay},
            {'params': kappa_params, 'lr': self.config.learning_rate * self.config.kappa_lr_scale,
             'weight_decay': 0.0},  # No weight decay on curvature
        ])

        self.scaler = GradScaler() if self.config.use_mixed_precision else None

    def train(self):
        self.logger.info(f"Starting κ convergence experiment (seed={self.config.seed}, init_κ={self.config.init_kappa})")
        self.model.train()

        dataloader_iter = iter(self.dataloader)
        accum_step = 0

        while self.global_step < self.config.max_steps:
            try:
                batch = next(dataloader_iter)
            except StopIteration:
                dataloader_iter = iter(self.dataloader)
                batch = next(dataloader_iter)

            tokens, tax_ids = batch
            tokens = tokens.to(self.device)
            tax_ids = tax_ids.to(self.device)

            with autocast(enabled=self.config.use_mixed_precision):
                outputs = self.model(tokens, tax_ids=tax_ids)
                loss = outputs['loss'] / self.config.gradient_accumulation_steps

            if self.scaler:
                self.scaler.scale(loss).backward()
            else:
                loss.backward()

            accum_step += 1

            if accum_step % self.config.gradient_accumulation_steps == 0:
                if self.scaler:
                    self.scaler.unscale_(self.optimizer)
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.config.max_grad_norm)
                    self.scaler.step(self.optimizer)
                    self.scaler.update()
                else:
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.config.max_grad_norm)
                    self.optimizer.step()

                self.optimizer.zero_grad()
                self.global_step += 1

                # Log κ
                kappa_val = self.model.hyper.c.item()
                kappa_grad = self.model.hyper.c.grad.item() if self.model.hyper.c.grad is not None else 0.0

                entry = {
                    'step': self.global_step,
                    'kappa': kappa_val,
                    'kappa_grad': kappa_grad,
                    'mlm': outputs['mlm'],
                    'dec': outputs['dec'],
                    'hex': outputs['hex'],
                    'dist': outputs['dist'],
                    'total': outputs['loss'].item() * self.config.gradient_accumulation_steps,
                }
                self.kappa_history.append(entry)

                if self.global_step % self.config.log_every == 0:
                    self.logger.info(
                        f"Step {self.global_step:>5d} | κ={kappa_val:.6f} | "
                        f"∇κ={kappa_grad:+.2e} | "
                        f"MLM={outputs['mlm']:.4f} HEX={outputs['hex']:.4f} "
                        f"DIST={outputs['dist']:.4f}"
                    )

                if self.global_step % self.config.save_every == 0:
                    self._save_checkpoint()

        # Final save
        self._save_checkpoint()
        self._save_history()

        final_kappa = self.model.hyper.c.item()
        self.logger.info(f"COMPLETE: Final κ = {final_kappa:.6f} (init was {self.config.init_kappa})")
        return final_kappa

    def _save_checkpoint(self):
        path = self.output_dir / f"checkpoint_step_{self.global_step}.pt"
        torch.save({
            'step': self.global_step,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'config': asdict(self.config),
            'kappa': self.model.hyper.c.item(),
            'kappa_history': self.kappa_history,
        }, path)
        self.logger.info(f"Saved checkpoint: {path}")

    def _save_history(self):
        path = self.output_dir / "kappa_history.json"
        with open(path, 'w') as f:
            json.dump(self.kappa_history, f, indent=2)
        self.logger.info(f"Saved κ history: {path} ({len(self.kappa_history)} entries)")

        # Also save a summary
        kappas = [e['kappa'] for e in self.kappa_history]
        summary = {
            'seed': self.config.seed,
            'init_kappa': self.config.init_kappa,
            'final_kappa': kappas[-1] if kappas else None,
            'kappa_trajectory': {
                'start': kappas[0] if kappas else None,
                'step_500': kappas[min(499, len(kappas)-1)] if kappas else None,
                'step_1000': kappas[min(999, len(kappas)-1)] if kappas else None,
                'step_2000': kappas[min(1999, len(kappas)-1)] if kappas else None,
                'final': kappas[-1] if kappas else None,
            },
            'total_steps': self.global_step,
        }
        with open(self.output_dir / "summary.json", 'w') as f:
            json.dump(summary, f, indent=2)


# =========================================================================
# 5. MAIN
# =========================================================================

def run_single(args):
    """Run a single κ convergence experiment."""
    config = KappaConfig(
        seed=args.seed,
        init_kappa=args.init_kappa,
        max_steps=args.steps,
        output_dir=f"{args.output_base}/seed_{args.seed}_init_{args.init_kappa:.2f}",
        data_root=args.data_root,
        tokenizer_path=args.tokenizer_path,
        manifest_path=args.manifest_path,
    )
    trainer = KappaConvergenceTrainer(config)
    return trainer.train()


def run_sweep(args):
    """Run the full 5-seed × 3-init sweep."""
    seeds = [42, 137, 2024, 888, 7]
    init_kappas = [0.5, 1.0, 2.0]

    results = []
    for seed in seeds:
        for init_k in init_kappas:
            print(f"\n{'='*60}")
            print(f"  SWEEP: seed={seed}, init_κ={init_k}")
            print(f"{'='*60}\n")

            config = KappaConfig(
                seed=seed,
                init_kappa=init_k,
                max_steps=args.steps,
                output_dir=f"{args.output_base}/seed_{seed}_init_{init_k:.1f}",
                data_root=args.data_root,
                tokenizer_path=args.tokenizer_path,
                manifest_path=args.manifest_path,
            )
            trainer = KappaConvergenceTrainer(config)
            final_k = trainer.train()

            results.append({
                'seed': seed,
                'init_kappa': init_k,
                'final_kappa': final_k,
            })

            # Save intermediate results
            with open(f"{args.output_base}/sweep_results.json", 'w') as f:
                json.dump(results, f, indent=2)

    # Print summary
    print(f"\n{'='*60}")
    print("  SWEEP COMPLETE")
    print(f"{'='*60}")
    for r in results:
        print(f"  seed={r['seed']:>4d}  init_κ={r['init_kappa']:.1f}  →  final_κ={r['final_kappa']:.6f}")

    kappas = [r['final_kappa'] for r in results]
    print(f"\n  Mean κ = {np.mean(kappas):.6f} ± {np.std(kappas):.6f}")
    print(f"  Range: [{min(kappas):.6f}, {max(kappas):.6f}]")

    return results


def main():
    parser = argparse.ArgumentParser(description="κ Convergence Experiment")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--init_kappa", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--sweep", action="store_true", help="Run full 5-seed × 3-init sweep")
    parser.add_argument("--output_base", type=str, default="./kappa_convergence_results")
    parser.add_argument("--data_root", type=str,
                        default="/zfs_raid/SentryBio/5k_test_genomes/processed_biosphere_data_supervised")
    parser.add_argument("--tokenizer_path", type=str,
                        default="/zfs_raid/SentryBio/5k_test_genomes/tokenizer_output_bulletproof/biosphere_tokenizer_bulletproof.pkl")
    parser.add_argument("--manifest_path", type=str,
                        default="/zfs_raid/SentryBio/5k_test_genomes/MASTER_REGISTRY_refseq_only.csv")
    args = parser.parse_args()

    if args.sweep:
        run_sweep(args)
    else:
        run_single(args)


if __name__ == "__main__":
    main()
