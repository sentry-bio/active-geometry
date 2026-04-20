# Validation Suite

This directory contains all empirical validation experiments for the
Active Geometry theory.

## Overview

| Validation | Method | Key Result |
|------------|--------|------------|
| **Genomic** (§4.1) | 5-seed training at κ=1.0 fixed | Procrustes r = 0.94 ± 0.02 |
| **Phylogenetic** (§4.2) | Direct H² tree embedding | κ = 3.0 (fungi), 12.7 (archaea), 16.4 (bacteria), all n ≈ 2 |
| **Viral** (§5) | 15 RNA virus families | Pearson r = 0.996 with state-equation prediction |
| **Protein** (§6) | 15 Pfam families | κ = 3.80 ± 0.60 vs. predicted 3.90 (2.6% agreement) |
| **Falsification** (§5.3) | 3 negative controls | Euclidean κ=0.00, synthetic MRE 1.08%, destroyed-structure r<0.3 |

## Directory Structure

```
validation/
├── genomic/           # DNA training experiments
│   ├── scripts/       # Analysis code
│   └── results/       # YAML result files
├── viral/             # RNA virus validation
│   ├── scripts/       # Sweep and analysis code
│   └── results/       # YAML result files
└── phylogenetic/      # Tree-based validation
    ├── scripts/       # κ estimation pipeline
    ├── trees/         # Newick/Nexus files
    └── results/       # YAML result files
```

## Running Validations

### Genomic (requires GPU, 6-10 hours)

```bash
# From repository root
make train-all-seeds
make analyze-convergence
```

### Viral (requires GPU, 2-4 hours)

```bash
make viral-validation
```

### Phylogenetic (CPU only, ~30 minutes)

```bash
make tree-validation
```

## Results Format

All results are stored as structured YAML for machine readability:

```yaml
experiment:
  name: "..."
  date_measured: "..."

summary:
  mean_procrustes_r: 0.9482
  std_procrustes_r: 0.0163
  # or for tree-embedding results:
  # kappa: 3.0
  # n_backsolved: 2.00
```

These files are cross-checked against `constants.yaml` by the
`Validate Scientific Consistency` CI workflow.
