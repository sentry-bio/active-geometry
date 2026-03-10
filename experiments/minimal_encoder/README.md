# Minimal Encoder Validation (SI Section 12)

Architecture-independent validation of the coordinate system using a
40,482-parameter minimal encoder:

- Single convolutional layer + 3-layer MLP
- No classification heads, no ODE flow
- Maps directly to H^2(kappa = 5/4) with curvature fixed analytically
- Training signal: quartet consistency (NCBI taxonomy) + radial ordering (genome size)

## Results

Five independent seeds yield mean Procrustes residual of 0.020 across 268
organisms. The sole undetermined degree of freedom is a global SO(2) rotation.

## Usage

```bash
python train.py \
    --manifest path/to/manifest.csv \
    --output-dir ./results \
    --device cuda \
    --n-seeds 5
```

## Files

- `train.py` — Training script (self-contained)

Coordinate tables and training logs are generated in the output directory.
