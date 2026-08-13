# Formal proofs in Lean 4

This directory machine-checks the algebraic spine of Active Geometry. The
authoritative mathematical treatment, including the metric packing argument
and its hypotheses, is [`../MATHEMATICAL_SPINE.md`](../MATHEMATICAL_SPINE.md).

## Mathematical hierarchy

The principal coordinate-free statement is the addressability bound

\[
\beta\le c\,h_{\mathrm{vol}},
\]

where:

- \(\beta\) is retained-information growth in nats per generative step;
- \(c\) converts generative steps to radial distance;
- \(h_{\mathrm{vol}}\) is host volume entropy in nats per radial distance.

For an isotropic hyperbolic host,
\(h_{\mathrm{vol}}=(n-1)\sqrt\kappa\), so

\[
\kappa\ge
\left(\frac{\beta}{c(n-1)}\right)^2.
\]

Capacity saturation gives equality. If
\(\beta=h_{\mathrm{eff}}\ln2\) and
\(\bar\kappa:=c^2\kappa\), the normalized equality is

\[
\bar\kappa=
\left(\frac{h_{\mathrm{eff}}\ln2}{n-1}\right)^2.
\]

The familiar formula without \(c\) is raw curvature only in the
process-time gauge \(c=1\).

## Files

```text
ActiveGeometry/
├── Addressability.lean    # scale-aware bound and normalized equality
└── KappaCurvature.lean    # normalized state-equation algebra and ceilings
```

### `Addressability.lean`

Formalized results include:

| Declaration | Meaning |
|---|---|
| `addressability_forces_positive_entropy` | \(\beta>0\), \(c>0\), and \(\beta\le c h_{\rm vol}\) imply \(h_{\rm vol}>0\) |
| `efficiency_le_one` | \(\eta=\beta/(c h_{\rm vol})\le1\) |
| `curvature_at_least_floor` | isotropic addressability gives a curvature lower bound |
| `saturated_curvature_eq_floor` | saturation fixes raw curvature once \(c,n\) are fixed |
| `floor_saturates_capacity` | the curvature floor realizes equality |
| `normalized_floor_eq_ideal` | multiplying by \(c^2\) removes radial-unit dependence |
| `process_time_gauge` | \(c=1\) recovers the familiar formula |
| `normalized_curvature_scale_invariant` | \(c^2\kappa\) is invariant under radial rescaling |

### `KappaCurvature.lean`

This file retains the normalized algebra:

- positivity and uniqueness of the equality value;
- monotonicity in transmitted rate and ambient dimension;
- alphabet-capacity ceilings;
- non-negative rate-matching potentials and their unique zero.

Its declaration `κ h n` now denotes ideal **normalized** curvature
\(\bar\kappa\), not unit-independent raw sectional curvature.

## What Lean does not establish

The formalization intentionally does not claim to machine-check:

1. the metric packing theorem that derives
   \(\beta\le c h_{\mathrm{vol}}\);
2. equivalence of packing and volume entropy under bounded geometry;
3. the space-form classification or the hyperbolic volume formula;
4. the Buneman/Gromov tree-classification theorems;
5. Sarkar's low-distortion embedding theorem;
6. Fisher–Rao curvature computations;
7. a physical dynamics toward capacity saturation;
8. empirical membership of any biological or linguistic system.

These are respectively paper proofs, classical cited results, open modeling
choices, or empirical questions.

The functions `U` and `V` in `KappaCurvature.lean` are positive-definite
mismatch functions. Their non-negativity and unique-zero theorems do not prove
Lyapunov stability without an explicit evolution law and a proof that the
function decreases along its trajectories.

## Tree dimension

The four-point condition classifies exact tree metrics. Minimal smooth ambient
dimension \(n=2\) comes from embeddability: a genuinely branching tree cannot
live faithfully in a connected one-dimensional Riemannian manifold, while
finite trees admit arbitrarily low-distortion embeddings in
\(\mathbb H^2\).

The Lean theorem that normalized curvature decreases with \(n\) is algebraic
monotonicity. It is not a proof that an objective selects \(n=2\).

## Build

```bash
cd theory/lean
lake build
```

Requires Lean 4 and Mathlib. The checked files contain no `sorry` declarations.
