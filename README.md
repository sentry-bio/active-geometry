# Active Geometry

**The Intrinsic Hyperbolic Curvature of Evolution**

[![Lean 4](https://img.shields.io/badge/Lean-4-blue)](theory/lean/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)

> Evolution organizes genomic sequences on a hyperbolic manifold with universal curvature **κ = 1.247 ± 0.003**, matching the theoretical prediction **κ = 1.23** with 1.7% agreement and zero free parameters.

---

## Overview

This repository provides complete, reproducible materials for two companion papers:

| Paper | Focus | Key Result |
|-------|-------|------------|
| **Empirical** | Measurement & Validation | κ = 1.247 ± 0.003 from 5,627 genomes |
| **Theory** | First-Principles Derivation | κ = (h ln 2)² ≈ 1.23 with n=2 |

Both papers share:
- Formal proofs in [Lean 4](theory/lean/)
- Canonical constants in [`constants.yaml`](constants.yaml)
- The BiosphereCodec model in [`model/`](model/)

---

## Quick Start

### One-Command Verification (Docker)

```bash
# Build and run all verifications
docker build -t active-geometry .
docker run --rm active-geometry
```

This executes:
1. Lean 4 formal proofs (machine-checked)
2. All 6 validation notebooks (numerical)
3. Constants consistency check

### Manual Verification

#### 1. Verify Lean Proofs (5 minutes)

```bash
cd theory/lean
lake build
```

#### 2. Reproduce κ Measurement (6-10 hours GPU)

```bash
# Train 5 independent models
make train-all-seeds

# Analyze coordinate convergence
make analyze-convergence
```

### 3. Validate on RNA Viruses (2-4 hours GPU)

```bash
# Run 15-dataset validation sweep
make viral-validation
```

### 4. Validate on Phylogenetic Trees (30 minutes CPU)

```bash
# κ estimation from tree geometry
make tree-validation
```

### 5. Generate All Figures

```bash
make figures
```

---

## Repository Structure

```
active-geometry/
│
├── manifest.yaml               # Machine-readable claims & verification map
├── constants.yaml              # Single source of truth for all values
├── Dockerfile                  # One-command reproducibility
├── run_all_verifications.sh    # Verification orchestrator
│
├── theory/
│   └── lean/                   # Formal proofs (Lean 4, machine-checked)
│       ├── ActiveGeometry/
│       │   └── KappaCurvature.lean  # 577 lines of proofs
│       └── lakefile.lean
│
├── model/                      # BiosphereCodec implementation
│   ├── biosphere_codec.py      # Core encoder-decoder
│   ├── training.py             # Training pipeline
│   └── hyperbolic.py           # Poincaré geometry utilities
│
├── validation/
│   ├── notebooks/              # 6 canonical verification notebooks
│   │   ├── 01_neural_convergence.ipynb
│   │   ├── 02_theory_verification.ipynb
│   │   ├── 03_viral_validation.ipynb
│   │   ├── 04_null_simulations.ipynb
│   │   ├── 05_topology.ipynb
│   │   └── 06_ablation.ipynb
│   ├── genomic/                # DNA training (5,627 genomes)
│   ├── viral/                  # RNA validation (15 viruses)
│   └── phylogenetic/           # Tree κ estimation
│
├── supplementary/
│   └── wolfram/                # CAS verification (optional, commercial)
│
├── data/
│   └── manifests/              # Data acquisition specs
│
└── figures/
    └── fig*.py                 # Reproducible figure generation
```

---

## The Geometric State Equation

The central theoretical result is a **parameter-free prediction** of evolutionary curvature:

$$\kappa = \left(\frac{h \ln 2}{n-1}\right)^2$$

Where:
- **κ** = Gaussian curvature of the Poincaré ball
- **h** = Entropy rate of genomic sequences (~1.6 bits/nucleotide)
- **n** = Intrinsic dimensionality (empirically n = 2)

Substituting measured values:
$$\kappa = \left(\frac{1.61 \times 0.693}{2-1}\right)^2 = 1.23$$

This matches the empirical measurement **κ = 1.247 ± 0.003** within 1.7%.

---

## Key Results

### Empirical (5-seed training)

| Seed | κ | Convergence |
|------|---|-------------|
| 0 | 1.245 | Reference |
| 42 | 1.248 | r = 0.984 |
| 123 | 1.247 | r = 0.981 |
| 456 | 1.249 | r = 0.979 |
| 789 | 1.246 | r = 0.983 |
| **Mean** | **1.247 ± 0.003** | **CV = 0.24%** |

### Viral Validation (15 datasets)

- κ correlates with phylogenetic depth: **ρ = 0.84, p < 0.001**
- Young viruses (SARS-CoV-2, Influenza): κ ≈ 1.32-1.35
- Ancient viruses (HCV, DENV): κ ≈ 1.35-1.55
- Validates substrate-independence (RNA vs DNA)

### Formal Verification

All theoretical claims verified in **Lean 4**:
- State equation derivation
- Dimensional analysis
- Uniqueness of n = 2 solution

---

## Installation

### Requirements

- Python 3.10+
- PyTorch 2.0+
- CUDA 11.8+ (for GPU training)
- Lean 4 (for proof verification)

### Setup

```bash
# Clone repository
git clone https://github.com/sentry-bio/active-geometry.git
cd active-geometry

# Create environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
make test
```

### Optional: Lean 4

```bash
# Install elan (Lean version manager)
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh

# Build proofs
cd theory/lean && lake build
```

---

## Reproducibility

### Data Access

Raw genomic data is not included due to size. To replicate:

```bash
# Download 1,540 public NCBI genomes (Tier 1)
make fetch-data

# Or use provided manifests to fetch specific accessions
python data/scripts/fetch_from_manifest.py --manifest data/manifests/public_refseq.tsv
```

### Checkpoints

Pre-trained checkpoints for verification:

```bash
# Download reference checkpoint (seed 42)
make download-checkpoint
```

### Validation

```bash
# Run all validation checks
make validate-all

# This runs:
# 1. Lean proof compilation
# 2. Python unit tests
# 3. Results consistency check (vs constants.yaml)
# 4. Figure regeneration
```

---

## Citation

If you use this work, please cite both papers:

```bibtex
@article{fenn2026empirical,
  title={Evolution as Active Geometry: A Universal Curvature Constant},
  author={Fenn, Rohit and Fenn, Amit},
  journal={},
  year={2026}
}

@article{fenn2026theory,
  title={A Geometric State Equation for Evolutionary Dynamics},
  author={Fenn, Rohit and Fenn, Amit},
  journal={},
  year={2026}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

- Repository: [github.com/sentry-bio/active-geometry](https://github.com/sentry-bio/active-geometry)
- Issues: [GitHub Issues](https://github.com/sentry-bio/active-geometry/issues)
