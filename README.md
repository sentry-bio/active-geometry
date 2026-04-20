# Active Geometry

**The Geometric State Equation of the Tree of Life**

[![Lean 4](https://img.shields.io/badge/Lean-4-blue)](theory/lean/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)
[![Preprint](https://img.shields.io/badge/bioRxiv-2026.03.09.710612-red.svg)](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2)

> Evolution is two-dimensional. Across every system tested — from
> decade-old viral outbreaks to 3.8-billion-year cellular lineages,
> across DNA, RNA, and protein alphabets — the intrinsic embedding
> dimension is **n = 2.00 ± 0.05**. Curvature is scale-dependent and
> predicted from the entropy rate of the generating code by the state
> equation **κ = (h · ln 2 / (n−1))²**, with zero adjustable parameters.
> The universal invariant is the dimension, not the curvature.

---

## What This Repository Contains

Reproducible materials for Fenn & Fenn (2026), [*Evolution as Active
Geometry: The Geometric State Equation of the Tree of Life*](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2).

- **Lean 4 proofs** of the state equation's mathematical consequences
  (existence, uniqueness, monotonicity, maximization at n=2, Lyapunov
  stability). Machine-checked, zero sorries.
- **BiosphereCodec reference encoder**: a minimal open-source Poincaré-ball
  encoder that reproduces the paper's Procrustes coordinate convergence
  on 5,550 genomes. (The production-scale encoder, Biosphere Atlas v8.4
  with 125K species coverage, is described in a forthcoming applications
  paper; this is the reproducibility kernel.)
- **Six validation notebooks** covering neural convergence, theory
  verification, viral families, null simulations, topology sweep,
  and ablation.
- **Encoder-free methods**: §4.2 κ-sweep tree embedding has no neural
  network — any reviewer can reproduce κ = 3–16 on GTDB/Li2021 trees
  with standard MDS + L-BFGS.

---

## The Two Findings

The paper reports two results of different character:

**1. Topological invariant — robust across everything tested.**

Back-solving `n = 1 + h·ln 2/√κ` from independently measured (h, κ)
pairs returns n = 2 at every scale, for every alphabet, across 10⁶-fold
variation in timescale and 10⁴-fold variation in mutation rate. This
is the deepest result.

| System | Substrate | κ | n |
|---|---|---|---|
| Multi-domain cellular life | DNA | 1.28–1.34 | 2.01 |
| Viral families (15) | RNA/DNA | 1.32–1.55 | 2.00 ± 0.05 |
| Fungi (Li 2021) | DNA | 3.0 ± 0.1 | 2.00 |
| Archaea (GTDB r220) | DNA | 12.7 ± 0.6 | 1.99 |
| Bacteria (GTDB r220) | DNA | 16.4 ± 0.5 | 1.99 |
| Protein families (15 Pfam) | Amino acid | 3.80 ± 0.60 | 2.03 ± 0.10 |

**2. Metric law — scale-dependent, predicted by entropy.**

The state equation `κ = (h·ln 2/(n−1))²` applied at the scale-appropriate
h predicts κ across a 13-fold range, with no parameters fit. Cross-alphabet
validation: h_protein = 2.85 bits predicts κ_protein = 3.90, measured
3.80 ± 0.60 — a 3.1× jump from the DNA value, confirmed at 2.6% agreement.

Viral curvature-entropy correlation: **Pearson r = 0.996**, explaining
99.3% of variance with zero free parameters.

---

## Quick Start

### Encoder-free verification (30 minutes, CPU only)

The single strongest empirical result in the paper — §4.2 domain-level
tree embedding — uses no neural network:

```bash
make tree-validation
```

This reproduces Table 3: fungi κ = 3.0, archaea κ = 12.7, bacteria
κ = 16.4, all with n = 2.

### Lean proofs (5 minutes)

```bash
cd theory/lean && lake build
```

Compiles all theorems including `kappa_n2`, `kappa_max_at_n2`,
`kappa_bounded_by_alphabet`, `lyapunov_zero_iff`, and the rest. Zero
sorries. See [theory/lean/README.md](theory/lean/README.md) for the full
theorem inventory.

### Docker (one-command, reproducible)

```bash
docker build -t active-geometry .
docker run --rm active-geometry
```

Runs Lean proofs + all 6 validation notebooks + constants consistency check.

### Five-seed Procrustes convergence (6–10 hours GPU)

```bash
make train-all-seeds
make analyze-convergence
```

Trains five independent BiosphereCodec instances with κ fixed at 1.0,
measures Procrustes alignment across all 10 seed pairs. Expected result:
mean r = 0.94 ± 0.02 (paper §4.1, SI §3).

---

## Repository Structure

```
active-geometry/
├── manifest.yaml              # Machine-readable claims & verification map
├── constants.yaml             # Single source of truth (scale-stratified κ, h, n)
├── Dockerfile                 # One-command reproducibility
├── run_all_verifications.sh   # Verification orchestrator
│
├── theory/lean/               # Lean 4 formalization
│   └── ActiveGeometry/KappaCurvature.lean
│
├── model/                     # BiosphereCodec reference encoder
│   ├── biosphere_codec.py
│   ├── training.py
│   └── hyperbolic.py
│
├── validation/
│   ├── notebooks/             # Six canonical notebooks
│   ├── genomic/               # 5,550-genome training + Procrustes
│   ├── viral/                 # 15-family curvature-entropy sweep
│   └── phylogenetic/          # Domain-level tree embedding + telescopes
│
├── supplementary/wolfram/     # Optional CAS cross-verification
├── data/manifests/            # Public genome accessions (RefSeq subset)
└── figures/                   # Reproducible figure generation
```

---

## What Is and Isn't Claimed

**Is claimed.** The state equation is a parameter-free law that holds
across biological systems tested. Dimension is the universal invariant.
The neural encoder's coordinate system is intrinsic to the data up to
global SO(2) rotation.

**Is not claimed.** Curvature does not emerge from training by gradient
descent on sequence data — the paper explicitly states (§7.5) that
curvature is set as a design parameter from theory, not discovered.
In the reference five-seed experiment (§4.1, SI §3), κ is fixed at
1.0 and what converges is coordinate geometry, measured by Procrustes
alignment (r = 0.94 ± 0.02). The curvature value κ ≈ 1.25 comes from
theory (state equation at h = 1.61) and from post-hoc telescope sweeps
on frozen embeddings (κ ≈ 1.28–1.34).

**Falsification criteria.** §7.5 names four conditions under which the
framework is refuted. Three negative controls have been run — Euclidean
null (κ = 0), synthetic recovery (1.08% MRE), destroyed structure
(Procrustes r < 0.3). All pass. The fourth — independent replication
on different architectures — is partially supported by Pearce et al.
2025 (Evo 2 DNA LM with sparse autoencoder analysis independently finds
phylogenetic geometry as a curved manifold).

---

## Citation

```bibtex
@article{fenn2026evolution,
  title={Evolution as Active Geometry: The Geometric State Equation of the Tree of Life},
  author={Fenn, Rohit and Fenn, Amit},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.03.09.710612}
}
```

---

## License

MIT. See [LICENSE](LICENSE).

## Contact

- Preprint: [bioRxiv 2026.03.09.710612](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2)
- Correspondence: research@sentry.bio
- Issues: [GitHub Issues](https://github.com/sentry-bio/active-geometry/issues)
