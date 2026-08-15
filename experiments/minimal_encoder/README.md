# Minimal Encoder Validation (SI Section 12)

Architecture-independent **Layer IIa instrument**: a 40,482-parameter
encoder that inhabits polar \(\mathbb H^2\), rather than discovering a
host by fitting curvature.

- Single convolutional layer + 3-layer MLP
- No classification heads, no ODE flow
- Maps directly to \(\mathbb H^2\) with \(\kappa\) **frozen** (the
  reference run uses \(5/4\); that number is a design choice, not a
  genetic-code theorem — InfoNCE temperature is degenerate with
  curvature, so \(\kappa\) cannot be learned)
- Two training axes, kept apart on purpose:
  - quartet consistency from NCBI taxonomy (angular / genealogical splits)
  - radial ordering from genome size (a depth proxy, **not** accumulated
    information and not a clock)

This is the embeddability floor \(n=2\): one radial coordinate for
process depth, one angular coordinate for divergence. Path trees remain
one-dimensional. The quartet loss is Theorem 6.1 used as a training
signal — classification, not curvature calibration.

## What a positive result supports

Five independent seeds yield mean Procrustes residual of 0.020 across 268
organisms. The sole undetermined degree of freedom is a global
\(\mathrm{SO}(2)\) rotation.

That is seed-stable **host-class** evidence: the same polar chart, up to
gauge. It does not certify saturation, an absolute \(\kappa\), or a
filled atlas of life. Genome-size radius stays advisory (E6). The two
axes are a better independence split than reading both from an inferred
tree, and still not a representation metric independent of taxonomy.

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
