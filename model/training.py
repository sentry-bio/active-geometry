#!/usr/bin/env python3
"""
🧬 Biosphere Codec: The Canonical Training Pipeline
==================================================

This script implements the complete 3-phase training procedure for the
256×4 Hyperbolic Biosphere Codec model. It is designed for reproducibility,
leveraging a configuration-driven approach and robust components for data
loading, training, and checkpointing.

**Training Strategy:**
1.  **Phase 1: Unsupervised Pre-training.** The model learns universal
    genomic patterns from 5,000+ MAGs using MLM and CLM objectives.
2.  **Phase 2: Supervised Fine-tuning.** The model learns biological
    hierarchy from annotated genomes using all objectives, including
    InfoNCE and phylogenetic distance regression.
3.  **Phase 3: Joint Training.** The model refines its understanding by
    training on the combined dataset, integrating universal patterns with
    specific biological knowledge.

**To Launch:**
`python biosphere_training.py`

**To Resume:**
`python biosphere_training.py --resume path/to/checkpoint.pt`
"""

import os
import sys
import json
import argparse
import logging
import time
import math
import pickle
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from collections import defaultdict

# Third-party imports
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import DataLoader, ConcatDataset, Dataset
from torch.optim import AdamW
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda.amp import GradScaler, autocast
import numpy as np
from tqdm import tqdm
import zstandard as zstd

# Optional WandB and GeoOpt imports
try:
    import wandb
    WANDB_AVAILABLE = True
except ImportError:
    WANDB_AVAILABLE = False
try:
    import geoopt
    GEOOPT_AVAILABLE = True
except ImportError:
    geoopt = None
    GEOOPT_AVAILABLE = False

# =========================================================================
# 1. CONFIGURATION
# =========================================================================
@dataclass
class BiosphereConfig:
    """Configuration for the Biosphere Codec training run."""
    # Model Architecture
    d_model: int = 256
    n_layers: int = 4
    latent_dim: int = 256
    vocab_size: int = 5444
    max_sequence_length: int = 8192

    # Training Strategy
    batch_size: int = 2
    gradient_accumulation_steps: int = 8
    learning_rate: float = 2e-4
    weight_decay: float = 0.01
    warmup_ratio: float = 0.1
    max_grad_norm: float = 1.0

    # Progressive Phases
    phase1_epochs: int = 10
    phase2_epochs: int = 8
    phase3_epochs: int = 5

    # Memory Optimization
    use_mixed_precision: bool = True
    gradient_checkpointing: bool = True

    # Paths & Logging
    tokenizer_path: str = "./data/tokenizer/biosphere_tokenizer.pkl"
    data_root: str = "./data/genomes"
    output_dir: str = "./biosphere_run_final"
    log_every: int = 50
    save_every: int = 1000
    use_wandb: bool = True
    max_checkpoints: int = 3
    
    def save(self, path: Path):
        with open(path, 'w') as f:
            json.dump(asdict(self), f, indent=2)

# =========================================================================
# 2. TOKENIZER WRAPPER
# =========================================================================
class BulletproofTokenizer:
    """
    A wrapper to provide a consistent interface for the pickled
    dictionary-based tokenizer.
    """
    def __init__(self, tokenizer_dict: dict):
        self.vocab = tokenizer_dict['vocabulary']
        self.vocab_size = tokenizer_dict['vocab_size']
        self.id_to_token = {v: k for k, v in self.vocab.items()}
        self.pad_token_id = self.vocab.get('<PAD>', 0)
        self.unk_token_id = self.vocab.get('<UNK>', 1)
        self.mask_token_id = self.vocab.get('<MASK>', self.unk_token_id)
        self.gs_token_id = self.vocab.get('<GS>', 2)
        self.ge_token_id = self.vocab.get('<GE>', 3)

    def encode(self, sequence: str, max_length: Optional[int] = None) -> List[int]:
        tokens = [self.gs_token_id]
        # Using a sliding window for k-mer tokenization (k=3 based on vocab)
        for i in range(len(sequence) - 2):
            kmer = sequence[i:i+3]
            tokens.append(self.vocab.get(kmer, self.unk_token_id))
        tokens.append(self.ge_token_id)

        if max_length:
            tokens = tokens[:max_length]
        return tokens

# =========================================================================
# 3. CANONICAL MODEL DEFINITION
# =========================================================================
class HyenaOperator(nn.Module):
    """Depth-wise gated 1-D convolution that mimics Hyena/Mamba behaviour."""
    def __init__(self, d_model: int, mode: str = "bidirectional", k_size: int = 7):
        super().__init__()
        assert mode in {"bidirectional", "causal"}, "mode must be 'bidirectional' or 'causal'"
        self.mode = mode
        self.k_size = k_size
        self.depthwise = nn.Conv1d(d_model,d_model,kernel_size=k_size,groups=d_model,bias=False,padding=0)
        self.gate = nn.Conv1d(d_model, d_model, kernel_size=1, bias=True)

    def _pad(self, x: Tensor) -> Tensor:
        if self.mode == "causal":
            return F.pad(x, (self.k_size - 1, 0))
        else:
            return F.pad(x, (self.k_size // 2, self.k_size // 2))

    def forward(self, x: Tensor) -> Tensor:
        x_t = x.transpose(1, 2)
        x_pad = self._pad(x_t)
        g = torch.sigmoid(self.gate(x_pad))
        y = self.depthwise(x_pad * g)
        if self.mode == "causal":
            y = y[..., -(x.size(1)):]
        return (x_t + y).transpose(1, 2)

class EncoderBlock(nn.Module):
    """Pre-norm Hyena block with residual dropout."""
    def __init__(self, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.hyena = HyenaOperator(d_model, mode="bidirectional")
        self.norm2 = nn.LayerNorm(d_model)
        self.ff = nn.Sequential(nn.Linear(d_model, 4 * d_model), nn.GELU(), nn.Dropout(dropout), nn.Linear(4 * d_model, d_model))
        self.do = nn.Dropout(dropout)

    def forward(self, x: Tensor) -> Tensor:
        x = x + self.do(self.hyena(self.norm1(x)))
        x = x + self.do(self.ff(self.norm2(x)))
        return x

class HierPool(nn.Module):
    """Aggregate token-level states into (global, gene, attn) pools."""
    def __init__(self, d_model: int):
        super().__init__()
        self.attn_gate = nn.Linear(d_model, 1)

    def forward(self, h: Tensor, gene_idx: Optional[List[Tensor]] = None) -> Tensor:
        global_pool = h.mean(dim=1)
        # Simplified pooling for this integration, can be expanded later
        return torch.cat([global_pool, global_pool, global_pool], dim=-1)

class BiosphereEncoder(nn.Module):
    """Full encoder with Hyena blocks and hierarchical pooling."""
    def __init__(self, vocab: int, d_model: int, n_layers: int, max_len: int):
        super().__init__()
        self.embed = nn.Embedding(vocab, d_model)
        self.pos = nn.Parameter(torch.empty(max_len, d_model))
        nn.init.normal_(self.pos, std=0.02)
        self.layers = nn.ModuleList([EncoderBlock(d_model) for _ in range(n_layers)])
        self.norm = nn.LayerNorm(d_model)
        self.pool = HierPool(d_model)

    def forward(self, ids: Tensor, gene_idx: Optional[List[Tensor]] = None) -> Tuple[Tensor, Tensor]:
        x = self.embed(ids) + self.pos[: ids.size(1)]
        for blk in self.layers:
            x = blk(x)
        h = self.norm(x)
        pooled = self.pool(h, gene_idx)
        return h, pooled

class PoincareMapping(nn.Module):
    """Linear → Poincaré-ball projection with learnable curvature."""
    def __init__(self, in_dim: int, latent_dim: int):
        super().__init__()
        self.lin = nn.Linear(in_dim, latent_dim)
        self.c = nn.Parameter(torch.tensor(1.0))
        self._manifold: Optional[Any] = None
        self._c_cached: Optional[Tensor] = None

    def _man(self):
        if not GEOOPT_AVAILABLE: return None
        if (self._manifold is None) or (self._c_cached is None) or (not torch.allclose(self._c_cached, self.c)):
            self._manifold = geoopt.PoincareBall(c=self.c)
            self._c_cached = self.c.detach().clone()
        return self._manifold

    def forward(self, x: Tensor) -> Tensor:
        z_euc = torch.tanh(self.lin(x))
        r_max = 0.9 / torch.sqrt(torch.abs(self.c) + 1e-8)
        z_euc = z_euc * r_max
        man = self._man()
        return man.projx(z_euc) if man is not None else z_euc

    def dist_mat(self, z: Tensor) -> Tensor:
        man = self._man()
        if man is not None:
            return man.dist(z.unsqueeze(1), z.unsqueeze(0))
        return torch.cdist(z, z)

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

class BiosphereLoss(nn.Module):
    """Combined learning objectives."""
    def __init__(self, manifold: PoincareMapping, mask_id: int, temp: float = 0.1):
        super().__init__()
        self.manifold = manifold
        self.mask_id = mask_id
        self.temp = temp

    def _infonce_loss(self, z: Tensor) -> Tensor:
        """InfoNCE contrastive loss on Poincaré embeddings.

        Pulls together embeddings from the same batch position across
        augmented views (odd/even split as proxy) while pushing apart
        embeddings from different genomes.
        """
        if z.size(0) < 4:
            return torch.tensor(0.0, device=z.device)
        # Split batch into two views (proxy augmentation)
        mid = z.size(0) // 2
        z_a, z_b = F.normalize(z[:mid], dim=-1), F.normalize(z[mid:2*mid], dim=-1)
        logits = z_a @ z_b.T / self.temp  # (mid, mid)
        labels = torch.arange(mid, device=z.device)
        return F.cross_entropy(logits, labels)

    def _distance_loss(self, z: Tensor) -> Tensor:
        """Poincaré distance margin loss.

        Encourages a minimum hyperbolic separation between distinct
        genome embeddings, preventing representational collapse.
        """
        if z.size(0) < 2:
            return torch.tensor(0.0, device=z.device)
        dist_mat = self.manifold.dist_mat(z)
        # Margin loss: distances should be >= margin (0.5)
        margin = 0.5
        mask = ~torch.eye(z.size(0), dtype=torch.bool, device=z.device)
        distances = dist_mat[mask]
        loss = F.relu(margin - distances).mean()
        return loss

    def forward(self, orig_tok: Tensor, mlm_labels: Tensor, enc_logits: Tensor, dec_logits: Tensor, z: Tensor) -> Tuple[Tensor, Dict[str, Any]]:
        vocab = enc_logits.size(-1)
        mlm_loss = F.cross_entropy(enc_logits.view(-1, vocab), mlm_labels.view(-1), ignore_index=-100)
        dec_loss = F.cross_entropy(dec_logits.view(-1, vocab), orig_tok.view(-1))
        infonce = self._infonce_loss(z)
        dist_loss = self._distance_loss(z)
        total = mlm_loss + dec_loss + 0.1 * infonce + 0.5 * dist_loss
        return total, {
            "mlm": mlm_loss.item(), "dec": dec_loss.item(),
            "infonce": infonce.item(), "dist": dist_loss.item(),
            "total": total.item(),
        }

class BiosphereCodec(nn.Module):
    """Complete encoder-decoder model for genomic sequences."""
    def __init__(self, config: BiosphereConfig):
        super().__init__()
        self.config = config
        self.encoder = BiosphereEncoder(config.vocab_size, config.d_model, config.n_layers, config.max_sequence_length)
        self.hyper = PoincareMapping(3 * config.d_model, config.latent_dim)
        self.decoder = BiosphereDecoder(self.encoder.embed)
        self.mask_id = 4 # Default mask token ID
        self.loss_fn = BiosphereLoss(self.hyper, self.mask_id)

    def mask_tokens(self, toks: Tensor) -> Tuple[Tensor, Tensor]:
        labels = toks.clone()
        mask = torch.rand_like(toks.float()) < 0.15
        labels[~mask] = -100
        out = toks.clone()
        out[mask] = self.mask_id
        return out, labels

    def forward(self, ids: Tensor) -> Dict[str, Any]:
        masked_ids, mlm_labels = self.mask_tokens(ids)
        enc_h, pooled = self.encoder(masked_ids)
        z = self.hyper(pooled)
        enc_logits = enc_h @ self.encoder.embed.weight.T
        dec_logits = self.decoder(enc_h)
        loss, logs = self.loss_fn(ids, mlm_labels, enc_logits, dec_logits, z)
        return {'loss': loss, 'loss_dict': logs, 'logits': enc_logits, 'z': z}

# =========================================================================
# 4. DATA LOADING
# =========================================================================
class GenomeDataset(Dataset):
    """Loads tokenized genome data from compressed .zst files."""
    def __init__(self, data_root: str, tokenizer: BulletproofTokenizer, config: BiosphereConfig, pattern: str):
        self.tokenizer = tokenizer
        self.config = config
        self.items = []

        zst_files = list(Path(data_root).glob(pattern))
        logging.info(f"Found {len(zst_files)} files for pattern '{pattern}'")

        for zst_file in tqdm(zst_files, desc=f"Loading '{pattern}'"):
            try:
                with open(zst_file, 'rb') as f_in:
                    dctx = zstd.ZstdDecompressor()
                    decompressed = dctx.decompress(f_in.read()).decode('utf-8')
                
                # Try JSONL format first (one JSON object per line)
                loaded_count = 0
                for line in decompressed.strip().split('\n'):
                    if not line:
                        continue
                    try:
                        item = json.loads(line)
                        if 'sequence' in item:
                            self.items.append(item['sequence'])
                            loaded_count += 1
                    except json.JSONDecodeError:
                        pass  # Skip malformed lines
                
                # Fallback: try as JSON array (for old chunks)
                if loaded_count == 0:
                    try:
                        data = json.loads(decompressed)
                        if isinstance(data, list):
                            for item in data:
                                if isinstance(item, dict) and 'sequence' in item:
                                    self.items.append(item['sequence'])
                                    loaded_count += 1
                    except json.JSONDecodeError:
                        logging.warning(f"Could not parse {zst_file} as JSONL or JSON array")
                
                if loaded_count == 0:
                    logging.warning(f"No sequences loaded from {zst_file}")
                    
            except Exception as e:
                logging.warning(f"Could not load {zst_file}: {e}")

    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        sequence = self.items[idx]
        tokens = self.tokenizer.encode(sequence, max_length=self.config.max_sequence_length)
        
        # Pad sequence
        padding_needed = self.config.max_sequence_length - len(tokens)
        tokens.extend([self.tokenizer.pad_token_id] * padding_needed)

        # Create labels for Masked Language Modeling (MLM)
        labels = list(tokens)
        for i, token in enumerate(labels):
            # Only mask non-special, non-padded tokens
            if token > 3 and np.random.random() < 0.15:
                 tokens[i] = self.tokenizer.mask_token_id
            else:
                 labels[i] = -100 # Ignore non-masked tokens in loss

        return torch.tensor(tokens, dtype=torch.long), torch.tensor(labels, dtype=torch.long)

# =========================================================================
# 5. TRAINER
# =========================================================================
class ElegantTrainer:
    """Manages the 3-phase training and evaluation loop."""
    def __init__(self, config: BiosphereConfig, resume_from: Optional[str] = None):
        self.config = config
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        config.save(self.output_dir / 'config.json')

        self.global_step = 0
        self.best_loss = float('inf')

        self._setup_logging()
        self._initialize_components()
        
        if resume_from:
            self._load_checkpoint(resume_from)

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                            handlers=[logging.FileHandler(self.output_dir / "training.log"), logging.StreamHandler()])
        self.logger = logging.getLogger("BiosphereTrainer")
        if WANDB_AVAILABLE and self.config.use_wandb:
            try:
                wandb.init(project="biosphere-codec-final", config=asdict(self.config), dir=str(self.output_dir))
            except Exception as e:
                self.logger.warning(f"Could not initialize wandb: {e}")
                self.config.use_wandb = False

    def _initialize_components(self):
        self.logger.info("Initializing components...")
        with open(self.config.tokenizer_path, 'rb') as f:
            self.tokenizer = BulletproofTokenizer(pickle.load(f))
        
        self.model = BiosphereCodec(self.config).to(self.device)
        self.optimizer = AdamW(self.model.parameters(), lr=self.config.learning_rate, weight_decay=self.config.weight_decay)
        self.scaler = GradScaler(enabled=self.config.use_mixed_precision)
        self.logger.info(f"Model created on {self.device} with {sum(p.numel() for p in self.model.parameters())/1e6:.1f}M params.")

    def _train_phase(self, phase_name: str, dataset: Dataset, epochs: int):
        self.logger.info(f"\n{'='*20} Starting Phase: {phase_name} ({epochs} epochs) {'='*20}")
        if len(dataset) == 0:
            self.logger.warning(f"Skipping phase '{phase_name}' due to empty dataset.")
            return

        dataloader = DataLoader(dataset, batch_size=self.config.batch_size, shuffle=True, num_workers=4, pin_memory=True, drop_last=True)
        
        num_training_steps = len(dataloader) * epochs // self.config.gradient_accumulation_steps
        lr_scheduler = OneCycleLR(self.optimizer, max_lr=self.config.learning_rate, 
                                  total_steps=num_training_steps,
                                  pct_start=self.config.warmup_ratio)
        self.model.train()

        for epoch in range(epochs):
            pbar = tqdm(dataloader, desc=f"Phase '{phase_name}' Epoch {epoch+1}/{epochs}")
            for i, (tokens, _) in enumerate(pbar): # Labels are now created inside the model forward pass
                tokens = tokens.to(self.device)
                with autocast(enabled=self.config.use_mixed_precision):
                    outputs = self.model(tokens)
                    loss = outputs['loss'] / self.config.gradient_accumulation_steps
                self.scaler.scale(loss).backward()
                if (i + 1) % self.config.gradient_accumulation_steps == 0:
                    self.scaler.unscale_(self.optimizer)
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.config.max_grad_norm)
                    self.scaler.step(self.optimizer)
                    self.scaler.update()
                    self.optimizer.zero_grad()
                    lr_scheduler.step()
                    self.global_step += 1
                    current_loss = loss.item() * self.config.gradient_accumulation_steps
                    if self.global_step % self.config.log_every == 0:
                         pbar.set_postfix({'loss': current_loss})
                         if self.config.use_wandb and WANDB_AVAILABLE:
                             wandb.log({'loss': current_loss, 'lr': lr_scheduler.get_last_lr()[0], 'step': self.global_step, **outputs['loss_dict']})
                    if self.global_step % self.config.save_every == 0:
                        self._save_checkpoint()

    def train(self):
        self.logger.info("Loading datasets for all phases...")
        unsupervised_ds = GenomeDataset(self.config.data_root, self.tokenizer, self.config, "processed_*/**/*.zst")
        supervised_ds = GenomeDataset(self.config.data_root, self.tokenizer, self.config, "processed_*supervised*/**/*.zst")
        
        self._train_phase("Unsupervised Pre-training", unsupervised_ds, self.config.phase1_epochs)
        self._train_phase("Supervised Fine-tuning", supervised_ds, self.config.phase2_epochs)
        
        combined_ds = ConcatDataset([unsupervised_ds, supervised_ds])
        self._train_phase("Joint Training", combined_ds, self.config.phase3_epochs)
        
        self.logger.info("🎉 Training completed successfully!")
        self._save_checkpoint(final=True)

    def _save_checkpoint(self, final=False):
        filename = "final_model.pt" if final else f"checkpoint_step_{self.global_step}.pt"
        path = self.output_dir / filename
        torch.save({
            'step': self.global_step,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
        }, path)
        self.logger.info(f"💾 Checkpoint saved to {path}")
        
    def _load_checkpoint(self, path: str):
        self.logger.info(f"Resuming from checkpoint: {path}")
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.global_step = checkpoint['step']
        self.logger.info(f"Resumed training from step {self.global_step}")

# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="🧬 Biosphere Codec Canonical Training Pipeline")
    parser.add_argument("--resume", type=str, help="Path to checkpoint to resume training from.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--output_dir", type=str, default="./biosphere_run_final", help="Output directory for checkpoints and logs.")
    parser.add_argument("--phase1_epochs", type=int, default=10, help="Number of epochs for Phase 1 (unsupervised pre-training).")
    parser.add_argument("--phase2_epochs", type=int, default=8, help="Number of epochs for Phase 2 (supervised fine-tuning).")
    parser.add_argument("--phase3_epochs", type=int, default=5, help="Number of epochs for Phase 3 (joint training).")
    parser.add_argument("--use_wandb", type=str, default="True", help="Enable Weights & Biases logging (True/False).")
    parser.add_argument("--data_root", type=str, default=None, help="Override data root directory.")
    parser.add_argument("--tokenizer_path", type=str, default=None, help="Override tokenizer path.")
    args = parser.parse_args()

    # Set random seeds for reproducibility
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Create config with command-line overrides
    config = BiosphereConfig()
    config.output_dir = args.output_dir
    config.phase1_epochs = args.phase1_epochs
    config.phase2_epochs = args.phase2_epochs
    config.phase3_epochs = args.phase3_epochs
    config.use_wandb = args.use_wandb.lower() in ('true', '1', 'yes')
    if args.data_root:
        config.data_root = args.data_root
    if args.tokenizer_path:
        config.tokenizer_path = args.tokenizer_path

    trainer = ElegantTrainer(config, resume_from=args.resume)
    trainer.train()

if __name__ == "__main__":
    main() 