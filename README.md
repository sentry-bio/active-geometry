# Active Geometry

**The Addressability Limit**

[![Lean 4](https://img.shields.io/badge/Lean-4-blue)](theory/lean/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)
[![Preprint](https://img.shields.io/badge/bioRxiv-2026.03.09.710612-red.svg)](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2)

> Remembering while creating costs room. Retained-history growth obeys the
> coordinate-free limit **β ≤ c·h_pack**, and block capacity reaches that
> ceiling exactly. Distinguishable endpoints need not preserve their
> genealogy; in a real hyperbolic host, however, genealogy has zero
> exponential-order tax.

The one-page dependency graph is
[`theory/THROUGHLINE.md`](theory/THROUGHLINE.md): **limit → ladder → quartet
classifier → polar chart**. It joins the mathematical balloon (addresses
without genealogy) to the surviving canonical-coordinate-system construction
(\(\mathbb H^2\) as depth plus divergence), without requiring nature to
saturate.

The formal kernel is in
[`theory/ADDRESSABILITY_KERNEL.md`](theory/ADDRESSABILITY_KERNEL.md); the full
proof, units, scope, and falsification criteria are in
[`theory/MATHEMATICAL_SPINE.md`](theory/MATHEMATICAL_SPINE.md). The inequality
is the general limit. The state equation is its capacity-saturating, isotropic
equality case; raw curvature follows only after a radial gauge is fixed.

---

## What This Repository Contains

Reproducible materials for Fenn & Fenn (2026), [*Evolution as Active
Geometry: The Geometric State Equation of the Tree of Life*](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2).

- **Lean 4 proofs** of the finite-depth metric packing count, the
  convergent-rate addressability theorem, scale-aware curvature floor,
  normalized scale-invariance, and conditional equality. Machine-checked,
  zero sorries. The Lean files do not claim the full limsup extension,
  space-form classification, or a physical dynamics.
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

Those \(\kappa\) and \(n\) rows below are the preprint's reported figure.
They are not certified measurements, and they are not the throughline
([`theory/THROUGHLINE.md`](theory/THROUGHLINE.md)).

The paper reports two results of different character:

**1. Topological classification — tree-like descent has minimal ambient n=2.**

The four-point condition classifies exact tree metrics, and finite trees admit
low-distortion embeddings in `H²`. Back-solving
`n = 1 + h·ln 2/(c√κ)` from independently measured quantities is a consistency
check, not an independent derivation of `n = 2`. Empirically, the reported
systems return values near two across 10⁶-fold variation in timescale and
10⁴-fold variation in mutation rate.

| System | Substrate | κ | n |
|---|---|---|---|
| Multi-domain cellular life | DNA | 1.28–1.34 | 2.01 |
| Viral families (15) | RNA/DNA | 1.32–1.55 | 2.00 ± 0.05 |
| Fungi (Li 2021) | DNA | 3.0 ± 0.1 | 2.00 |
| Archaea (GTDB r220) | DNA | 12.7 ± 0.6 | 1.99 |
| Bacteria (GTDB r220) | DNA | 16.4 ± 0.5 | 1.99 |
| Protein families (15 Pfam) | Amino acid | 3.80 ± 0.60 | 2.03 ± 0.10 |

**2. Metric limit — host capacity is bounded below by retained-information growth.**

The general statement is `β ≤ c·h_vol`. For an isotropic hyperbolic host this
becomes the curvature floor
`κ ≥ (h·ln 2/(c(n−1)))²`; an economical, capacity-saturating representation
achieves equality. Cross-alphabet measurements reported by the paper include
`h_protein = 2.85` bits predicting normalized `κ̄_protein = 3.90`, with measured
raw `κ_protein = 3.80 ± 0.60` under the paper's radial convention.

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

Compiles the metric packing theorem and scale-aware addressability algebra,
including `faithful_representation_addressable`,
`curvature_at_least_floor`, `normalized_state_equation`,
`normalized_curvature_scale_invariant`, and the existing alphabet-capacity
bounds. Zero sorries. See
[theory/lean/README.md](theory/lean/README.md) for scope and inventory.

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
├── theory/
│   ├── ADDRESSABILITY_KERNEL.md # compact theorem/assumption boundary
│   ├── MATHEMATICAL_SPINE.md  # Definitions, packing proof, equality conditions
│   └── lean/                  # Lean 4 algebraic formalization
│       └── ActiveGeometry/
│           ├── Packing.lean
│           ├── Addressability.lean
│           └── KappaCurvature.lean
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

**Is claimed mathematically.** Retained exponential novelty obeys
`β ≤ c·h_pack`; for Riemannian hosts with uniform fixed-ball volume bounds,
`h_pack = h_vol`. Polynomial-growth hosts cannot support positive retained-
information growth at finite radial rate. Exact tree metrics satisfy the
four-point condition; genuinely branching examples in the chosen smooth
hyperbolic embedding class have minimal ambient dimension two.

**Is tested empirically.** Biological tree metrics prefer hyperbolic
representations; reported independently estimated rates and curvatures are
compared with the capacity-saturating equality case. The neural encoder's
coordinate system is reproducible across seeds up to global rotation.

**Is not claimed.** Curvature does not emerge from training by gradient
descent on sequence data — the paper explicitly states (§7.5) that
curvature is set as a design parameter from theory, not discovered.
In the reference five-seed experiment (§4.1, SI §3), κ is fixed at
1.0 and what converges is coordinate geometry, measured by Procrustes
alignment (r = 0.94 ± 0.02). A positive mismatch potential alone does not
establish a physical dynamics or Lyapunov attraction toward equality.

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
