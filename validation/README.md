# Validation Suite

This directory contains all empirical validation experiments for the
Active Geometry theory.

## Overview

| Validation | Method | Key Result |
|------------|--------|------------|
| **Genomic** | 5-seed training | κ = 1.247 ± 0.003 |
| **Viral** | 15 RNA virus sweeps | ρ = 0.84 with phylogenetic depth |
| **Phylogenetic** | Tree embedding | κ = 1.245 ± 0.015 (model-free) |

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
  date: "..."

summary:
  kappa_mean: 1.247
  kappa_std: 0.003
```

These files are validated against `constants.yaml` by CI.
