#!/usr/bin/env python3
"""
κ Convergence Experiment — Canonical Reproduction v2
====================================================

Fixed version of run_kappa_convergence.py, informed by the full model
evolution record (BIOSPHERE_MODEL_EVOLUTION.md, E6 v5→v10b trajectory).

Three bugs in v1 prevented κ from receiving gradient:
  1. AMP (autocast float16) killed gradients through torch.acosh
  2. Phylum-level positives at batch_size=8 → ~0 positive pairs per batch
  3. Margin DIST loss saturated to 0 early, contributing nothing to κ

The fix applies the E6 v10b lesson: κ converges to 5/4 if and only if the
loss makes κ load-bearing for pairwise metric structure. InfoNCE does this
because distance RATIOS (not just orderings) depend on curvature. CE on
prototypes does not (argmin is curvature-invariant).

Design:
  - Full BiosphereCodec architecture (Hyena encoder + Poincaré head)
  - MLM + DEC for encoder gradients (they do NOT inform κ — confirmed)
  - Genus-level InfoNCE on Poincaré distances (THE κ-informing loss)
  - κ reparameterized via softplus for clean gradients
  - Separate κ optimizer: Adam lr=1e-2, no weight decay
  - NO mixed precision (float32 throughout — acosh needs it)
  - Batch size 32 with genus-balanced sampling (ensures positive pairs)

Expected result: κ converges to ~1.25 regardless of seed or init,
with ∇κ showing actual gradient signal (not 0.00).

Usage:
  python run_kappa_convergence_v2.py --sweep
  python run_kappa_convergence_v2.py --seed 42 --init_kappa 0.5 --steps 5000
"""

import os
import sys
import json
import argparse
import logging
import math
import time
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import DataLoader, Dataset, Sampler, WeightedRandomSampler
from torch.optim import Adam, AdamW
import numpy as np


# =========================================================================
# 1. CONFIGURATION
# =========================================================================

@dataclass
class KappaConfig:
    """Configuration for κ convergence experiment."""
    # Architecture (matches BPE tokenizer: vocab=4096, tokens are 512-d npy)
    d_model: int = 256
    n_layers: int = 4
    latent_dim: int = 128
    vocab_size: int = 4096     # BPE tokenizer (v10_1_tokenized)
    max_sequence_length: int = 512  # pre-tokenized to 512 tokens per genome

    # Curvature
    init_kappa: float = 1.0

    # Training
    batch_size: int = 32       # 32 fits 6GB VRAM; grouped sampler guarantees pos pairs
    members_per_family: int = 2  # 2 = max diversity (16 fam × 2), 4 = more pairs (8 × 4)
    sampler_mode: str = 'stochastic'  # 'stochastic' = WeightedRandom (Nexus), 'grouped' = deterministic
    gradient_accumulation_steps: int = 2  # effective batch = 128
    learning_rate: float = 2e-4
    kappa_lr: float = 1e-2     # Separate, higher lr for κ (E6 v10b: 1e-3 with AdamW)
    weight_decay: float = 0.01
    max_grad_norm: float = 1.0
    max_steps: int = 5000

    # Loss weights
    hex_weight: float = 0.5    # v1 was 0.1 — InfoNCE needs to be dominant for κ signal
    infonce_temp: float = 0.1  # τ=0.1 amplifies distance ratio sensitivity 10×
    mlm_prob: float = 0.15

    # Data — uses v5 manifest with pre-tokenized .npy files
    manifest_path: str = "/zfs_raid/SentryBio/working/ultimate_training_manifest_v5_enriched_tokenized_full.csv"
    max_genomes: int = 5000    # sample ~5k from 44k manifest
    min_family_count: int = 5  # Minimum genomes per family for InfoNCE positives

    # Logging
    seed: int = 42
    output_dir: str = "./kappa_convergence_v2"
    log_every: int = 10
    save_every: int = 1000


# =========================================================================
# 2. MODEL COMPONENTS (from BiosphereCodec, minimal extraction)
# =========================================================================

class HyenaOperator(nn.Module):
    """Hyena convolution operator."""
    def __init__(self, d_model: int, mode: str = "encoder"):
        super().__init__()
        self.depthwise = nn.Conv1d(d_model, d_model, 7, padding=3, groups=d_model)
        self.gate = nn.Linear(d_model, d_model)

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
        pooled = torch.cat([mean_pool, max_pool, attn_pool], dim=-1)

        return h, pooled


class PoincareMapping(nn.Module):
    """Linear → Poincaré ball projection with LEARNABLE curvature.

    v2 fix: softplus reparameterization for clean gradients.
    c = softplus(raw_c), so c is always positive and ∂c/∂raw_c = sigmoid(raw_c).
    No geoopt needed.
    """
    def __init__(self, in_dim: int, latent_dim: int = 128, init_c: float = 1.0):
        super().__init__()
        self.lin = nn.Linear(in_dim, latent_dim)
        # Softplus reparameterization: c = softplus(raw_c)
        # raw_c = softplus_inverse(init_c) = log(exp(init_c) - 1)
        if init_c > 0.01:
            raw_init = math.log(math.exp(init_c) - 1.0)
        else:
            raw_init = -4.0
        self.raw_c = nn.Parameter(torch.tensor(raw_init))

    @property
    def c(self) -> Tensor:
        """Actual curvature, always positive via softplus."""
        return F.softplus(self.raw_c)

    def _project_to_ball(self, z: Tensor) -> Tensor:
        """Project points inside the Poincaré ball of curvature c."""
        c = self.c
        r_max = 0.9 / torch.sqrt(c)
        norm = torch.norm(z, dim=-1, keepdim=True).clamp(min=1e-8)
        max_norm = r_max - 1e-5
        factor = torch.where(norm > max_norm, max_norm / norm, torch.ones_like(norm))
        return z * factor

    def forward(self, x: Tensor) -> Tensor:
        z = torch.tanh(self.lin(x))
        c = self.c
        r_max = 0.9 / torch.sqrt(c)
        z = z * r_max
        return self._project_to_ball(z)

    def poincare_dist(self, u: Tensor, v: Tensor) -> Tensor:
        """Poincaré ball distance with gradients through c.

        d(u,v) = (1/√c) · arccosh(1 + 2c‖u-v‖² / ((1-c‖u‖²)(1-c‖v‖²)))

        This is the critical path: every distance depends on c,
        so InfoNCE logits = -d/τ create a non-trivial curvature preference.
        """
        c = self.c
        eps = 1e-7

        sqdist = ((u - v) ** 2).sum(-1)
        sqnorm_u = (u ** 2).sum(-1)
        sqnorm_v = (v ** 2).sum(-1)

        denom = (1 - c * sqnorm_u) * (1 - c * sqnorm_v)
        denom = torch.clamp(denom, min=eps)

        x = 1 + 2 * c * sqdist / denom
        x = torch.clamp(x, min=1.0 + eps)

        dist = torch.acosh(x) / torch.sqrt(c + eps)
        return dist

    def dist_mat(self, z: Tensor) -> Tensor:
        """Full pairwise distance matrix [B, B]."""
        B = z.shape[0]
        u = z.unsqueeze(1).expand(B, B, -1)
        v = z.unsqueeze(0).expand(B, B, -1)
        return self.poincare_dist(u, v)


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
    """Full BiosphereCodec with learnable κ and genus-level InfoNCE.

    v2 changes from v1:
      - softplus κ reparameterization (clean gradients)
      - genus-level InfoNCE (not phylum — denser positive pairs)
      - no DIST margin loss (saturates to 0, replaced by InfoNCE)
      - no mixed precision (float32 for acosh stability)
    """

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

    def forward(self, ids: Tensor, family_ids: Optional[Tensor] = None) -> Dict[str, Any]:
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

        # --- MLM + DEC losses (encoder gradients only, NOT κ-informing) ---
        mlm_loss = F.cross_entropy(
            enc_logits.view(-1, vocab), labels.view(-1), ignore_index=-100
        )
        dec_loss = F.cross_entropy(dec_logits.view(-1, vocab), ids.view(-1))

        # --- HEX: Family-level InfoNCE on Poincaré distances ---
        # THIS is the κ-informing loss. Every distance passes through
        # poincare_dist(u, v, c=softplus(raw_c)), so ∂L/∂raw_c ≠ 0.
        #
        # v2 fix: use family (not phylum) for denser positive pairs.
        # With ~200 families having ≥5 members and family-balanced sampling,
        # a batch of 64 yields ~15-40 positive pairs per anchor.
        hex_loss = torch.tensor(0.0, device=ids.device)
        n_pos_pairs = 0

        if family_ids is not None and B >= 4:
            dist = self.hyper.dist_mat(z)          # [B, B]
            logits = -dist / self.config.infonce_temp  # similarity

            diag = torch.eye(B, dtype=torch.bool, device=ids.device)
            logits = logits.masked_fill(diag, -1e9)

            # Positive pairs: same family
            pos_mask = family_ids.unsqueeze(1) == family_ids.unsqueeze(0)
            pos_mask = pos_mask & ~diag

            if pos_mask.any():
                n_pos_pairs = pos_mask.sum().item()
                # Standard supervised InfoNCE:
                # For each anchor with positives, compute:
                #   -log(Σ exp(sim_pos)) + log(Σ exp(sim_all))
                has_pos = pos_mask.any(dim=1)

                if has_pos.sum() >= 2:
                    log_denom = torch.logsumexp(logits[has_pos], dim=1)  # [V]

                    # Mean positive log-prob per anchor
                    losses = []
                    for i, row_idx in enumerate(torch.where(has_pos)[0]):
                        pos_sims = logits[row_idx][pos_mask[row_idx]]
                        # log-mean-exp of positive similarities
                        log_num = torch.logsumexp(pos_sims, dim=0) - math.log(pos_sims.shape[0])
                        losses.append(-(log_num - log_denom[i]))

                    hex_loss = torch.stack(losses).mean()

        # --- Total loss ---
        # MLM + DEC drive encoder representations (do not inform κ)
        # HEX drives κ (every distance is a function of curvature)
        total = mlm_loss + dec_loss + self.config.hex_weight * hex_loss

        # Get κ gradient for logging (will be populated after backward)
        kappa_val = self.hyper.c.item()
        raw_c_grad = self.hyper.raw_c.grad.item() if self.hyper.raw_c.grad is not None else 0.0

        # Radius diagnostics: how close are embeddings to the ball boundary?
        with torch.no_grad():
            c_val = self.hyper.c
            r_max = (0.9 / torch.sqrt(c_val)).item()
            radii = torch.norm(z, dim=-1)  # [B]
            r_mean = radii.mean().item()
            r_std = radii.std().item()
            r_max_actual = radii.max().item()
            # Fraction of radius used (1.0 = at boundary)
            r_frac = r_mean / r_max if r_max > 0 else 0.0
            # Mean pairwise distance for scale context
            if family_ids is not None and B >= 4:
                dist_mean = dist.mean().item()
                dist_std = dist.std().item()
            else:
                dist_mean = 0.0
                dist_std = 0.0

        return {
            'loss': total,
            'mlm': mlm_loss.item(),
            'dec': dec_loss.item(),
            'hex': hex_loss.item(),
            'kappa': kappa_val,
            'raw_c': self.hyper.raw_c.item(),
            'kappa_grad': raw_c_grad,
            'n_pos_pairs': n_pos_pairs,
            'r_mean': r_mean,
            'r_std': r_std,
            'r_max': r_max_actual,
            'r_boundary': r_max,
            'r_frac': r_frac,
            'dist_mean': dist_mean,
            'dist_std': dist_std,
        }


# =========================================================================
# 3. DATASET
# =========================================================================

class GenomeDataset(Dataset):
    """Loads pre-tokenized .npy genomes with FAMILY labels for InfoNCE.

    Uses the v5 manifest (44k genomes, each pre-tokenized to 512 BPE tokens
    saved as .npy files). Samples ~max_genomes with family-balanced selection
    and filters to families with ≥min_family_count members.

    Family-level (not genus-level) because:
    - ~200-400 families vs ~1,100+ genera → denser positive pairs per batch
    - At batch_size=64: ~10-30 positive pairs vs ~2-4 with genus
    - Matches E6 v10b which used family-level InfoNCE
    """

    def __init__(self, manifest_path: str, max_genomes: int = 5000,
                 vocab_size: int = 4096, min_family_count: int = 5,
                 max_len: int = 512, seed: int = 42):
        import pandas as pd

        self.max_len = max_len
        self.vocab_size = vocab_size

        # Load manifest
        df = pd.read_csv(manifest_path, low_memory=False)
        df = df[df['tokenized_path'].notna() & df['family'].notna()]

        # Filter to files that exist (check a sample for speed)
        sample_check = df.sample(min(200, len(df)), random_state=seed)
        exist_rate = sample_check['tokenized_path'].apply(
            lambda p: os.path.exists(str(p))).mean()
        print(f"Manifest: {len(df)} genomes with tokenized_path + family "
              f"(~{exist_rate*100:.0f}% files exist)")
        if exist_rate < 0.5:
            # Full filter needed
            df = df[df['tokenized_path'].apply(lambda p: os.path.exists(str(p)))]
            print(f"After existence filter: {len(df)}")

        # Filter families with sufficient members for InfoNCE positives
        fc = df['family'].value_counts()
        valid_families = fc[fc >= min_family_count].index
        df = df[df['family'].isin(valid_families)]
        print(f"After family filter (>={min_family_count}): {len(df)} genomes, "
              f"{df['family'].nunique()} families")

        # Sample ~max_genomes with family-balanced selection
        if len(df) > max_genomes:
            n_families = df['family'].nunique()
            max_per_fam = max(min_family_count, max_genomes // n_families)
            sampled = []
            for fam, grp in df.groupby('family'):
                sampled.append(grp.sample(min(len(grp), max_per_fam), random_state=seed))
            df = pd.concat(sampled, ignore_index=True)
            if len(df) > max_genomes:
                df = df.sample(max_genomes, random_state=seed)
        print(f"Final sample: {len(df)} genomes, {df['family'].nunique()} families")

        self.paths = df['tokenized_path'].tolist()
        families = df['family'].tolist()
        self.family_to_id = {f: i for i, f in enumerate(sorted(set(families)))}
        self.family_labels = [self.family_to_id[f] for f in families]
        self.n_families = len(self.family_to_id)

        # Report distribution
        fam_counts = Counter(families)
        print(f"  {self.n_families} families, top 5 by count: "
              f"{fam_counts.most_common(5)}")
        # Expected positive pairs per batch of 64
        avg_fam_size = np.mean([n for n in fam_counts.values()])
        print(f"  Avg family size: {avg_fam_size:.1f} genomes")

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        tokens = np.load(self.paths[idx]).astype(np.int64)
        # Clamp to vocab range
        tokens = np.clip(tokens, 0, self.vocab_size - 1)
        # Pad or truncate to max_len
        if len(tokens) > self.max_len:
            tokens = tokens[:self.max_len]
        elif len(tokens) < self.max_len:
            tokens = np.pad(tokens, (0, self.max_len - len(tokens)))
        return (torch.from_numpy(tokens),
                torch.tensor(self.family_labels[idx], dtype=torch.long))

    def get_family_balanced_sampler(self) -> WeightedRandomSampler:
        """Family-balanced sampling ensures every batch has positive pairs."""
        label_counts = Counter(self.family_labels)
        weights = [1.0 / label_counts[label] for label in self.family_labels]
        return WeightedRandomSampler(weights, len(weights), replacement=True)

    def get_family_indices(self) -> Dict[int, List[int]]:
        """Return {family_id: [dataset indices]} for grouped sampling."""
        fam2idx: Dict[int, List[int]] = {}
        for idx, fam in enumerate(self.family_labels):
            fam2idx.setdefault(fam, []).append(idx)
        return fam2idx


class FamilyGroupedBatchSampler(Sampler):
    """Batch sampler that packs K families × M members per batch.

    Guarantees every batch contains positive pairs for InfoNCE.
    With batch_size=32, members_per_family=4: 8 families × 4 = 32 samples,
    giving 8 × C(4,2) = 48 positive pairs per batch.

    This is the fix for sparse-pairs on small GPUs where batch_size=64 won't fit.
    """

    def __init__(self, family_indices: Dict[int, List[int]],
                 batch_size: int = 32, members_per_family: int = 4,
                 seed: int = 42, drop_last: bool = True):
        self.batch_size = batch_size
        self.members_per_family = members_per_family
        self.families_per_batch = batch_size // members_per_family
        self.drop_last = drop_last
        self.rng = np.random.RandomState(seed)

        # Only keep families with enough members
        self.family_indices = {
            fam: idxs for fam, idxs in family_indices.items()
            if len(idxs) >= members_per_family
        }
        self.family_ids = list(self.family_indices.keys())

        # Estimate total batches per epoch
        total_samples = sum(len(v) for v in self.family_indices.values())
        self._len = total_samples // batch_size

    def __iter__(self):
        # Shuffle family order each epoch
        fam_order = self.family_ids.copy()
        self.rng.shuffle(fam_order)

        # Build an index pool per family (shuffled)
        pools = {}
        for fam in fam_order:
            idxs = self.family_indices[fam].copy()
            self.rng.shuffle(idxs)
            pools[fam] = idxs

        batches_yielded = 0
        fam_ptr = 0

        while fam_ptr + self.families_per_batch <= len(fam_order):
            batch = []
            for _ in range(self.families_per_batch):
                fam = fam_order[fam_ptr % len(fam_order)]
                fam_ptr += 1

                pool = pools[fam]
                # Sample members_per_family from this family (with wrap)
                for j in range(self.members_per_family):
                    if len(pool) == 0:
                        # Refill pool
                        pool = self.family_indices[fam].copy()
                        self.rng.shuffle(pool)
                        pools[fam] = pool
                    batch.append(pool.pop())

            yield batch
            batches_yielded += 1

        # If we exhausted families, wrap around for more batches
        # (ensures we don't under-train)

    def __len__(self):
        return self._len


# =========================================================================
# 4. TRAINER
# =========================================================================

class KappaConvergenceTrainer:
    """Trains BiosphereCodec with learnable κ and logs convergence.

    v2: No mixed precision, genus-level InfoNCE, softplus κ,
    separate κ optimizer, genus-balanced sampling.
    """

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
        self.logger = logging.getLogger("KappaConvergence_v2")

    def _init_model(self):
        self.model = BiosphereCodecKappa(self.config).to(self.device)
        n_params = sum(p.numel() for p in self.model.parameters())
        self.logger.info(f"Model: {n_params:,} parameters on {self.device}")
        self.logger.info(f"Initial κ = {self.model.hyper.c.item():.6f} "
                         f"(raw_c = {self.model.hyper.raw_c.item():.4f})")
        self.logger.info(f"raw_c requires_grad = {self.model.hyper.raw_c.requires_grad}")
        self.logger.info(f"NO mixed precision (float32 throughout)")

    def _init_data(self):
        dataset = GenomeDataset(
            manifest_path=self.config.manifest_path,
            max_genomes=self.config.max_genomes,
            vocab_size=self.config.vocab_size,
            min_family_count=self.config.min_family_count,
            max_len=self.config.max_sequence_length,
            seed=self.config.seed,
        )
        if self.config.sampler_mode == 'grouped':
            # Grouped batch sampler: deterministic K families × M members
            members_per_fam = self.config.members_per_family
            batch_sampler = FamilyGroupedBatchSampler(
                family_indices=dataset.get_family_indices(),
                batch_size=self.config.batch_size,
                members_per_family=members_per_fam,
                seed=self.config.seed,
            )
            self.dataloader = DataLoader(
                dataset, batch_sampler=batch_sampler,
                num_workers=2, pin_memory=True,
            )
            fams_per_batch = self.config.batch_size // members_per_fam
            pos_pairs = fams_per_batch * members_per_fam * (members_per_fam - 1)
            self.logger.info(f"Batch: {fams_per_batch} families × {members_per_fam} members = "
                             f"{self.config.batch_size} samples, ~{pos_pairs} guaranteed pos pairs")
        else:
            # Stochastic WeightedRandom sampler — the Nexus config that produced κ=1.25
            # Pos pairs emerge from random collisions, not structure.
            # This stochasticity may be essential for the correct κ basin.
            sampler = dataset.get_family_balanced_sampler()
            self.dataloader = DataLoader(
                dataset, batch_size=self.config.batch_size,
                sampler=sampler, num_workers=2, pin_memory=True, drop_last=True,
            )
            self.logger.info(f"Sampler: WeightedRandom (stochastic pos pairs)")

        self.logger.info(f"Dataset: {len(dataset)} genomes, {dataset.n_families} families")
        self.logger.info(f"Batch size: {self.config.batch_size}")

    def _init_optimizer(self):
        # CRITICAL: separate optimizer groups for κ vs encoder
        # κ gets higher lr (1e-2) with no weight decay
        # This follows E6 v10b which used AdamW lr=1e-3
        kappa_params = [self.model.hyper.raw_c]
        other_params = [p for n, p in self.model.named_parameters()
                        if p is not self.model.hyper.raw_c]

        self.optimizer = Adam([
            {'params': other_params, 'lr': self.config.learning_rate,
             'weight_decay': self.config.weight_decay},
            {'params': kappa_params, 'lr': self.config.kappa_lr,
             'weight_decay': 0.0},
        ])

        self.logger.info(f"Optimizer: Adam, encoder_lr={self.config.learning_rate}, "
                         f"κ_lr={self.config.kappa_lr}")

    def train(self):
        self.logger.info(f"Starting κ convergence v2 "
                         f"(seed={self.config.seed}, init_κ={self.config.init_kappa})")
        self.logger.info(f"HEX weight={self.config.hex_weight}, "
                         f"InfoNCE temp={self.config.infonce_temp}")
        self.model.train()

        dataloader_iter = iter(self.dataloader)
        accum_count = 0

        while self.global_step < self.config.max_steps:
            try:
                batch = next(dataloader_iter)
            except StopIteration:
                dataloader_iter = iter(self.dataloader)
                batch = next(dataloader_iter)

            tokens, family_ids = batch
            tokens = tokens.to(self.device)
            family_ids = family_ids.to(self.device)

            # NO autocast — float32 throughout for clean acosh gradients
            outputs = self.model(tokens, family_ids=family_ids)
            loss = outputs['loss'] / self.config.gradient_accumulation_steps

            loss.backward()
            accum_count += 1

            if accum_count % self.config.gradient_accumulation_steps == 0:
                # Clip and step — this is one optimizer step
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.config.max_grad_norm)

                # Log gradient BEFORE step (this is the accumulated gradient)
                raw_c_grad = (self.model.hyper.raw_c.grad.item()
                              if self.model.hyper.raw_c.grad is not None else 0.0)

                self.optimizer.step()
                self.optimizer.zero_grad()
                self.global_step += 1

                kappa_val = self.model.hyper.c.item()

                entry = {
                    'step': self.global_step,
                    'kappa': kappa_val,
                    'raw_c': self.model.hyper.raw_c.item(),
                    'kappa_grad': raw_c_grad,
                    'mlm': outputs['mlm'],
                    'dec': outputs['dec'],
                    'hex': outputs['hex'],
                    'n_pos_pairs': outputs['n_pos_pairs'],
                    'total': outputs['loss'].item() * self.config.gradient_accumulation_steps,
                }
                self.kappa_history.append(entry)

                if self.global_step % self.config.log_every == 0:
                    self.logger.info(
                        f"Step {self.global_step:>5d} | κ={kappa_val:.6f} | "
                        f"∇raw_c={raw_c_grad:+.4e} | "
                        f"MLM={outputs['mlm']:.4f} HEX={outputs['hex']:.4f} "
                        f"pos_pairs={outputs['n_pos_pairs']}"
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
            'raw_c': self.model.hyper.raw_c.item(),
            'kappa_history': self.kappa_history[-100:],  # last 100 for space
        }, path)
        self.logger.info(f"Saved checkpoint: {path}")

    def _save_history(self):
        path = self.output_dir / "kappa_history.json"
        with open(path, 'w') as f:
            json.dump(self.kappa_history, f, indent=2)
        self.logger.info(f"Saved κ history: {path} ({len(self.kappa_history)} entries)")

        # Summary
        kappas = [e['kappa'] for e in self.kappa_history]
        grads = [e['kappa_grad'] for e in self.kappa_history]
        summary = {
            'seed': self.config.seed,
            'init_kappa': self.config.init_kappa,
            'final_kappa': kappas[-1] if kappas else None,
            'kappa_trajectory': {
                'start': kappas[0] if kappas else None,
                'step_500': kappas[min(499, len(kappas)-1)] if kappas else None,
                'step_1000': kappas[min(999, len(kappas)-1)] if kappas else None,
                'step_2000': kappas[min(1999, len(kappas)-1)] if kappas else None,
                'step_3000': kappas[min(2999, len(kappas)-1)] if kappas else None,
                'final': kappas[-1] if kappas else None,
            },
            'gradient_stats': {
                'mean_abs_grad': float(np.mean(np.abs(grads))) if grads else None,
                'max_abs_grad': float(np.max(np.abs(grads))) if grads else None,
                'nonzero_frac': float(np.mean(np.array(grads) != 0.0)) if grads else None,
            },
            'total_steps': self.global_step,
            'hex_weight': self.config.hex_weight,
            'kappa_lr': self.config.kappa_lr,
            'theory': {
                'h': 1.61,
                'kappa_predicted': (1.61 * math.log(2)) ** 2,
            },
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
        manifest_path=args.manifest_path,
        max_genomes=args.max_genomes,
        latent_dim=args.latent_dim,
        min_family_count=args.min_family_count,
        batch_size=args.batch_size,
        members_per_family=args.members_per_family,
        sampler_mode=args.sampler_mode,
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
                manifest_path=args.manifest_path,
                max_genomes=args.max_genomes,
                latent_dim=args.latent_dim,
                min_family_count=args.min_family_count,
                batch_size=args.batch_size,
                members_per_family=args.members_per_family,
                sampler_mode=args.sampler_mode,
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
    mean_k = np.mean(kappas)
    std_k = np.std(kappas)
    cv = std_k / mean_k * 100 if mean_k > 0 else 999

    print(f"\n  Mean κ = {mean_k:.6f} ± {std_k:.6f} (CV = {cv:.1f}%)")

    # Theory comparison
    h = 1.61
    kappa_theory = (h * math.log(2)) ** 2
    agreement = abs(mean_k - kappa_theory) / kappa_theory * 100
    print(f"  Theory: κ = (h·ln2)² = ({h}×{math.log(2):.4f})² = {kappa_theory:.4f}")
    print(f"  Agreement: {agreement:.1f}%")

    if cv < 5:
        print(f"\n  ✓ CONVERGED: All runs agree within {cv:.1f}% CV")
    elif cv < 20:
        print(f"\n  ⚠ PARTIAL CONVERGENCE: CV = {cv:.1f}%")
        for ci in sorted(set(r['init_kappa'] for r in results)):
            subset = [r['final_kappa'] for r in results if r['init_kappa'] == ci]
            print(f"    init={ci:.1f}: {np.mean(subset):.6f} ± {np.std(subset):.6f}")
    else:
        print(f"\n  ✗ NO CONVERGENCE: CV = {cv:.1f}% — κ is init-dependent")

    return results


def main():
    parser = argparse.ArgumentParser(description="κ Convergence Experiment v2")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--init_kappa", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--sweep", action="store_true",
                        help="Run full 5-seed × 3-init sweep")
    parser.add_argument("--output_base", type=str,
                        default="./kappa_convergence_v2_results")
    parser.add_argument("--manifest_path", type=str,
                        default="/zfs_raid/SentryBio/working/ultimate_training_manifest_v5_enriched_tokenized_full.csv")
    parser.add_argument("--max_genomes", type=int, default=5000,
                        help="Max genomes to sample from manifest")
    parser.add_argument("--latent_dim", type=int, default=128,
                        help="Poincaré ball dimension (2 for H²)")
    parser.add_argument("--min_family_count", type=int, default=5,
                        help="Min genomes per family for InfoNCE positives")
    parser.add_argument("--batch_size", type=int, default=32,
                        help="Batch size (32 fits 6GB VRAM with grouped sampling)")
    parser.add_argument("--members_per_family", type=int, default=2,
                        help="Members per family in grouped sampler (2=max diversity, 4=more pairs)")
    parser.add_argument("--sampler_mode", type=str, default='stochastic',
                        choices=['stochastic', 'grouped'],
                        help="'stochastic' = WeightedRandom (Nexus config), 'grouped' = deterministic")
    args = parser.parse_args()

    if args.sweep:
        run_sweep(args)
    else:
        run_single(args)


if __name__ == "__main__":
    main()
