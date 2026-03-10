"""
Biosphere Codec
================
Minimal yet complete encoder-decoder for genomic sequences (~430 LOC).

Key features
-------------
1. Hyena-style state-space operator that runs in linear time and supports
   both bidirectional (encoder) and causal (decoder) modes.
2. Hierarchical pooling that aggregates token-level representations into
   gene-level and global sequence summaries, feeding a hyperbolic
   projection head.
3. Poincaré hyperbolic embedding with learnable curvature using *geoopt*
   when available (falls back to Euclidean distance if not installed).
4. Multi-task objective: masked-language modelling (MLM), causal LM,
   hierarchical InfoNCE (HEX) and optional patristic-distance regression.

This file is self-contained (no external Hyena, tokeniser, etc.) so that the
smoke-test at the bottom runs out-of-the-box. Replace the stub operator with
`hyena-ops` for better speed.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple, Any, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

# ---------------------------------------------------------------------------
# Optional hyperbolic library
# ---------------------------------------------------------------------------

try:
    import geoopt  # type: ignore
    _GEOOPT_AVAILABLE = True
except ImportError:  # pragma: no cover – runtime optional dependency
    geoopt = None  # type: ignore
    _GEOOPT_AVAILABLE = False


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 1. Hyena Operator (simplified, CPU-friendly)                             ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class HyenaOperator(nn.Module):
    """Depth-wise gated 1-D convolution that mimics Hyena/Mamba behaviour.

    *Bidirectional* variant uses symmetric *same* padding. The *causal*
    variant applies **left-only** padding to guarantee that output timestep *t*
    is a function of inputs `≤ t` only.
    """

    def __init__(self, d_model: int, mode: str = "bidirectional", k_size: int = 7):
        super().__init__()
        assert mode in {"bidirectional", "causal"}, "mode must be 'bidirectional' or 'causal'"

        self.mode = mode
        self.k_size = k_size

        # Depth-wise convolution (one kernel per channel)
        self.depthwise = nn.Conv1d(
            d_model,
            d_model,
            kernel_size=k_size,
            groups=d_model,
            bias=False,
            padding=0,  # we add padding manually in *forward*
        )

        # Channel-mixing *gate*
        self.gate = nn.Conv1d(d_model, d_model, kernel_size=1, bias=True)

    def _pad(self, x: Tensor) -> Tensor:  # x is [B, D, L]
        if self.mode == "causal":
            # left-pad with (k-1) zeros, keep right side intact
            # this ensures causal behavior - position t only sees positions ≤ t
            pad_amt = self.k_size - 1
            return F.pad(x, (pad_amt, 0))
        else:
            # symmetric *same* padding for bidirectional processing
            pad = self.k_size // 2
            return F.pad(x, (pad, pad))

    def forward(self, x: Tensor) -> Tensor:  # x: [B, L, D]
        x_t = x.transpose(1, 2)  # -> [B, D, L]

        x_pad = self._pad(x_t)
        # Important: Apply gate to padded signal so dimensions match
        g = torch.sigmoid(self.gate(x_pad))  # [B, D, L_pad]
        y = self.depthwise(x_pad * g)  # gated convolution

        if self.mode == "causal":
            # Strip the (k-1) left-most steps so output length == input length
            # This removes the extra padding we added on the left side
            y = y[..., -(x.size(1)):]

        return (x_t + y).transpose(1, 2)  # [B, L, D]


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 2. Encoder stack + Hierarchical pooling                                  ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class EncoderBlock(nn.Module):
    """Pre-norm Hyena block with residual dropout."""

    def __init__(self, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.hyena = HyenaOperator(d_model, mode="bidirectional")
        self.norm2 = nn.LayerNorm(d_model)

        self.ff = nn.Sequential(
            nn.Linear(d_model, 4 * d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(4 * d_model, d_model),
        )

        self.do = nn.Dropout(dropout)

    def forward(self, x: Tensor) -> Tensor:  # [B, L, D]
        x = x + self.do(self.hyena(self.norm1(x)))
        x = x + self.do(self.ff(self.norm2(x)))
        return x


class HierPool(nn.Module):
    """Aggregate token-level states into (global, gene, attn) pools."""

    def __init__(self, d_model: int):
        super().__init__()
        self.attn_gate = nn.Linear(d_model, 1)

    def forward(self, h: Tensor, gene_idx: Optional[List[Tensor]] = None) -> Tensor:
        """h: [B, L, D]; *gene_idx* is a list of 1-D index tensors (CPU or GPU)."""

        B, _, D = h.shape

        # ----- global mean -----
        global_pool = h.mean(dim=1)  # [B, D]

        # ----- gene-level mean -----
        if gene_idx is None:
            gene_pool = global_pool  # fall-back – still matches dimension
        else:
            gene_vecs: List[Tensor] = []
            for bi in range(B):
                # Handle case where gene_idx[bi] might be empty or None
                if gene_idx[bi] is None or gene_idx[bi].numel() == 0:
                    gene_vecs.append(h[bi].mean(0))
                    continue
                
                # Move indices to same device as h
                idx = gene_idx[bi].to(h.device)
                
                # Split sequence at gene boundaries
                splits = torch.tensor_split(h[bi], idx.tolist(), dim=0)
                
                # Calculate mean for each segment and then overall mean
                gene_means = torch.stack([seg.mean(0) for seg in splits if seg.numel() > 0])
                gene_vecs.append(gene_means.mean(0))
            
            gene_pool = torch.stack(gene_vecs, 0)  # [B, D]

        # ----- attention-weighted global -----
        w = torch.softmax(self.attn_gate(h), dim=1)  # [B, L, 1]
        attn_pool = (h * w).sum(dim=1)  # [B, D]

        return torch.cat([global_pool, gene_pool, attn_pool], dim=-1)  # [B, 3D]


class BiosphereEncoder(nn.Module):
    """Full encoder with Hyena blocks and hierarchical pooling."""
    
    def __init__(self, vocab: int, d_model: int = 512, n_layers: int = 6, max_len: int = 131_072):
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


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 3. Hyperbolic projection head                                           ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class PoincareMapping(nn.Module):
    """Linear → Poincaré-ball projection with learnable curvature using geoopt."""

    def __init__(self, in_dim: int, latent_dim: int = 128):
        super().__init__()
        self.lin = nn.Linear(in_dim, latent_dim)
        
        if _GEOOPT_AVAILABLE:
            # Initialize with κ=1.0 to test clean RNA curvature emergence
            initial_c = torch.tensor(1.0)
            self.manifold = geoopt.PoincareBall(c=initial_c, learnable=True)
        else:
            self.manifold = None
            # Fallback: simple learnable curvature parameter
            self.c = nn.Parameter(torch.tensor(1.56))

    def forward(self, x: Tensor) -> Tensor:  # [B, *]
        z_euc = torch.tanh(self.lin(x))  # bound to (-1, 1)
        
        if self.manifold is not None:
            # Use geoopt's intrinsic projection with learnable curvature
            # Scale to stay well inside the ball
            z_euc = z_euc * 0.9
            return self.manifold.projx(z_euc)
        else:
            # Euclidean fallback
            return z_euc

    def dist_mat(self, z: Tensor) -> Tensor:  # z: [B, d]
        if self.manifold is not None:
            return self.manifold.dist(z.unsqueeze(1), z.unsqueeze(0))  # [B, B]
        else:
            # Euclidean fallback
            return torch.cdist(z, z)
    
    @property
    def c(self) -> torch.Tensor:
        """Get current curvature parameter."""
        if self.manifold is not None and hasattr(self.manifold, 'c'):
            return self.manifold.c
        else:
            return getattr(self, '_c', torch.tensor(1.0))


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 4. Causal decoder                                                      ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class BiosphereDecoder(nn.Module):
    """Lightweight causal decoder with weight tying."""
    
    def __init__(self, shared_embed: nn.Embedding):
        super().__init__()
        d_model = shared_embed.embedding_dim

        self.hyena = HyenaOperator(d_model, mode="causal")
        self.norm = nn.LayerNorm(d_model)

        self.proj = nn.Linear(d_model, shared_embed.num_embeddings, bias=False)
        # Weight tying - same Parameter object
        self.proj.weight = shared_embed.weight

    def forward(self, h: Tensor) -> Tensor:  # [B, L, D]
        return self.proj(self.norm(self.hyena(h)))  # [B, L, V]


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 5. Multi-task loss                                                     ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class BiosphereLoss(nn.Module):
    """Combined learning objectives: MLM, CLM, InfoNCE, distance regression."""
    
    def __init__(self, manifold: PoincareMapping, mask_id: int, mlm_prob: float = 0.15, temp: float = 0.1, hex_weight: float = 0.1, dist_weight: float = 0.5):
        super().__init__()
        self.manifold = manifold
        self.mask_id = mask_id
        self.mlm_prob = mlm_prob
        self.temp = temp
        self.hex_weight = hex_weight
        self.dist_weight = dist_weight

    # ----------------------------- MLM masking helper

    def mask_tokens(self, toks: Tensor, vocab: int) -> Tuple[Tensor, Tensor]:
        """Create masked tokens and labels (-100 for unmasked positions)."""
        device = toks.device
        labels = toks.clone()

        mask = torch.rand_like(toks.float()) < self.mlm_prob
        labels[~mask] = -100  # ignore index for CE

        rand_tokens = torch.randint(0, vocab, toks.shape, device=device)
        mask_choice = torch.rand_like(toks.float())

        out = toks.clone()

        # 80 % → [MASK]
        out[mask & (mask_choice < 0.8)] = self.mask_id
        # 10 % → random token (0.8 ≤ mask_choice < 0.9)
        out[mask & ((mask_choice >= 0.8) & (mask_choice < 0.9))] = rand_tokens[mask & ((mask_choice >= 0.8) & (mask_choice < 0.9))]
        # 10 % → original token (already in place)

        return out, labels

    # ----------------------------- forward: compute losses

    def forward(
        self,
        orig_tok: Tensor,  # [B, L] original unmasked tokens
        mlm_labels: Tensor,  # [B, L] with -100 masking
        enc_logits: Tensor,  # [B, L, V]
        dec_logits: Tensor,  # [B, L, V]
        z: Tensor,  # [B, d]
        tax_ids: Optional[Tensor] = None,
        patristic: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Dict[str, Any]]:

        vocab = enc_logits.size(-1)

        # --- MLM cross-entropy (only on masked positions) ----------------
        mlm_loss = F.cross_entropy(
            enc_logits.view(-1, vocab), mlm_labels.view(-1), ignore_index=-100
        )

        # --- Causal LM on decoder ---------------------------------------
        dec_loss = F.cross_entropy(dec_logits.view(-1, vocab), orig_tok.view(-1))

        # --- Hierarchical InfoNCE ---------------------------------------
        hex_loss = torch.tensor(0.0, device=orig_tok.device)
        if tax_ids is not None:
            dist = self.manifold.dist_mat(z)  # [B,B]
            logits = -dist / self.temp  # similarity

            diag = torch.eye(len(z), dtype=torch.bool, device=z.device)
            logits = logits.masked_fill(diag, -1e9)  # exclude self

            levels = 1 if tax_ids.dim() == 1 else tax_ids.size(1)
            losses: List[Tensor] = []
            for lvl in range(levels):
                ids_lvl = tax_ids if tax_ids.dim() == 1 else tax_ids[:, lvl]
                pos_mask = ids_lvl.unsqueeze(1) == ids_lvl.unsqueeze(0)  # [B,B]
                pos_mask = pos_mask & ~diag
                
                if not pos_mask.any():
                    continue

                # InfoNCE per sample
                log_probs = logits - torch.logsumexp(logits, dim=1, keepdim=True)
                lvl_loss = -(log_probs[pos_mask].mean())
                losses.append(lvl_loss)

            if losses:
                hex_loss = torch.stack(losses).mean()

        # --- Distance regression ---------------------------------------
        dist_loss = torch.tensor(0.0, device=orig_tok.device)
        if patristic is not None:
            hyp = self.manifold.dist_mat(z)
            mask = torch.triu(torch.ones_like(hyp), diagonal=1).bool()
            dist_loss = F.mse_loss(hyp[mask], patristic[mask])

        total = mlm_loss + dec_loss + self.hex_weight * hex_loss + self.dist_weight * dist_loss

        return total, {
            "mlm": mlm_loss.item(),
            "dec": dec_loss.item(),
            "hex": hex_loss.item(),
            "dist": dist_loss.item(),
            "total": total.item(),
        }


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 6. Full model                                                          ║
# ╚═════════════════════════════════════════════════════════════════════════╝


class BiosphereCodec(nn.Module):
    """Complete encoder-decoder model for genomic sequences."""
    
    def __init__(
        self,
        vocab: int,
        d_model: int = 512,
        n_layers: int = 6,
        latent_dim: int = 128,
        mask_id: Optional[int] = None,
        hex_weight: float = 0.2,
        dist_weight: float = 0.5,
    ):
        super().__init__()

        # Core components
        self.encoder = BiosphereEncoder(vocab, d_model, n_layers)
        self.hyper = PoincareMapping(3 * d_model, latent_dim)
        self.decoder = BiosphereDecoder(self.encoder.embed)

        # Configure mask token ID (default to last token if not specified)
        self.mask_id = mask_id if mask_id is not None else vocab - 1
        self.loss_fn = BiosphereLoss(self.hyper, self.mask_id, hex_weight=hex_weight, dist_weight=dist_weight)

    # ---------------------------------------------------------------- encode

    def encode(self, ids: Tensor, gene_idx: Optional[List[Tensor]] = None) -> Tuple[Tensor, Tensor]:
        """Encode tokens to hidden states and latent hyperbolic codes."""
        h, pooled = self.encoder(ids, gene_idx)
        z = self.hyper(pooled)
        return h, z

    # ---------------------------------------------------------------- forward

    def forward(
        self,
        ids: Tensor,
        gene_idx: Optional[List[Tensor]] = None,
        tax_ids: Optional[Tensor] = None,
        patristic: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Dict[str, Any]]:
        """Full forward pass with loss computation."""

        # -------------- MLM masking -----------------------------------
        masked_ids, mlm_labels = self.loss_fn.mask_tokens(ids, self.encoder.embed.num_embeddings)

        # -------------- Encoder ---------------------------------------
        enc_h, z = self.encode(masked_ids, gene_idx)

        # Shared embedding projection for MLM logits (tied weights)
        enc_logits = enc_h @ self.encoder.embed.weight.T  # [B, L, V]

        # -------------- Decoder ---------------------------------------
        dec_logits = self.decoder(enc_h)

        # -------------- Loss ------------------------------------------
        # Important: Pass original tokens (ids) as first argument, not labels
        loss, logs = self.loss_fn(ids, mlm_labels, enc_logits, dec_logits, z, tax_ids, patristic)

        return loss, logs


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║ 7. Smoke-test                                                         ║
# ╚═════════════════════════════════════════════════════════════════════════╝


if __name__ == "__main__":
    torch.manual_seed(0)

    V = 4097  # 4096 DNA-BPE tokens + [MASK]
    model = BiosphereCodec(V, d_model=256, n_layers=2, latent_dim=32).eval()

    B, L = 3, 256
    ids = torch.randint(0, V - 1, (B, L))

    # Gene boundaries as list of tensors
    gene_idx: List[Tensor] = [torch.tensor([64, 128, 192]), torch.tensor([128]), torch.empty(0, dtype=torch.long)]

    tax = torch.randint(0, 10, (B, 2))  # genus + family (example)
    dist = torch.rand(B, B)

    loss, logs = model(ids, gene_idx, tax, dist)
    print("loss", loss.item())
    print(logs)
