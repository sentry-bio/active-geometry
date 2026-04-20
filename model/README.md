# BiosphereCodec Model

The core neural network for learning hyperbolic representations of genomic sequences.

## Architecture

```
BiosphereCodec (256×4)
├── Encoder
│   ├── Token Embedding (vocab → 256)
│   ├── Positional Encoding (learned, max 8192)
│   ├── Hyena Blocks × 4
│   │   ├── LayerNorm → HyenaOperator (bidirectional)
│   │   └── LayerNorm → FFN (256 → 1024 → 256)
│   └── Hierarchical Pooling → [global, gene, attention] → 768D
├── Poincaré Mapping
│   ├── Linear (768 → 128)
│   ├── Learnable curvature c (initialized to 1.0)
│   └── Ball projection (geoopt)
└── Decoder
    ├── HyenaOperator (causal)
    └── Weight-tied projection
```

## Key Features

1. **Hyena Operator**: Linear-time sequence mixing (no quadratic attention)
2. **Hierarchical Pooling**: Gene-aware aggregation
3. **Poincaré Curvature**: Set as a design parameter from the state equation
   κ = (h · ln 2)² (Paper I §7.5: curvature is degenerate with InfoNCE
   temperature in contrastive learning, so it cannot be discovered by
   gradient descent). The reference encoder trains at κ = 1.0 fixed; what
   is measured is coordinate-system stability via Procrustes alignment
   across seeds (Paper I §4.1, r = 0.94 ± 0.02).
4. **Multi-task Loss**: MLM + CLM + InfoNCE + Distance Regression

## Files

| File | Purpose |
|------|---------|
| `biosphere_codec.py` | Core model definition (~450 LOC) |
| `training.py` | 3-phase training pipeline |
| `hyperbolic.py` | Poincaré geometry utilities |

## Usage

```python
from model import BiosphereCodec

# Create model
model = BiosphereCodec(
    vocab=5444,
    d_model=256,
    n_layers=4,
    latent_dim=128
)

# Forward pass
tokens = torch.randint(0, 5444, (batch, seq_len))
loss, logs = model(tokens)

# Access the curvature parameter (fixed at training, not learned)
kappa = model.hyper.c.item()
print(f"Curvature: κ = {kappa:.3f}")
```

## Training

See [`training.py`](training.py) for the full pipeline:

```bash
python model/training.py --seed 42 --output_dir ./results
```

The training has 3 phases:
1. **Unsupervised**: MLM + CLM on raw sequences
2. **Supervised**: Add InfoNCE with taxonomic labels
3. **Joint**: Combined dataset refinement
