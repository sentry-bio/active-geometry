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

\(\mathbb H^2\) is the embeddability floor for genuinely branching trees
in the stated class. Interpreting radius as process depth and angle as
divergence is a modeling choice, not a consequence of the capacity
theorem. Path trees remain one-dimensional. The quartet loss is Theorem
6.1 used as a training signal — classification, not curvature
calibration.

## What a positive result supports

Five independent seeds yield mean Procrustes residual of 0.020 across 268
organisms. The leftover gauge is global \(O(2)\) (rotation and
reflection), unless an orientation convention is imposed.

That is seed-stable reproducibility **within the imposed model**:
\(\mathbb H^2\), frozen \(\kappa\), taxonomy quartets, and the radial
target are supplied by construction. It does not compare host classes
(occupancy of the growth-class \(\times\) tree-defect figure does; E9 is
a matched-packing illustration of Corollary 4.3). It does not certify saturation, an absolute \(\kappa\), or
a filled atlas of life. Genome-size radius stays advisory (E6). The two
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
