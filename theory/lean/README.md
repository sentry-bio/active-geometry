# Formal proofs in Lean 4

This directory machine-checks the metric packing count, its convergent-rate
addressability theorem, and the downstream algebra of the bound. The compact
dependency structure is
[`../ADDRESSABILITY_KERNEL.md`](../ADDRESSABILITY_KERNEL.md); the complete
metric packing argument and its hypotheses are in
[`../MATHEMATICAL_SPINE.md`](../MATHEMATICAL_SPINE.md).

The declaration mathematicians should cite is
`ActiveGeometry.addressability_limit`
(`Packing.addressability_limit`): in every proper metric host, a faithful
finite-resolution representation whose radii tend to infinity and whose
represented growth, radial rate, and packing growth converge satisfies
\(\beta\le c\,h_{\mathrm{pack}}\).

## Mathematical hierarchy

The theory has two layers (see
[`../ADDRESSABILITY_KERNEL.md`](../ADDRESSABILITY_KERNEL.md)). **Layer I** is
universal and curvature-free: the addressability bound, the block-capacity
identity, and the constrained-capacity ladder. **Layer II** is the curvature
realization: a space-form identification of \(h_{\mathrm{cap}}\) gives a
curvature floor; saturation makes that floor an equality (the face of the
feasible region). Lean formalizes the packing/block portion of Layer I and the
algebraic skeleton of Layer II. The growth-class identities in
`Measurability.lean` are instrument mathematics for the occupancy gate
(see [`../MEASURABILITY.md`](../MEASURABILITY.md)), not a third layer.

The weighted relational-capacity lower bound (Theorem 4.4, Skenderi) is a
paper proof, not Lean. Theorem 7.1 (Heintze isotropy / axiom A3) is a paper
sketch, not Lean.

The principal coordinate-free statement is the Layer I addressability bound

\[
\beta\le c\,h_{\mathrm{pack}},
\]

where:

- \(\beta\) is represented-history growth in nats per generative step;
- \(c\) converts generative steps to radial distance;
- \(h_{\mathrm{pack}}\) is host packing entropy in nats per radial distance.

For a host whose packing rate equals the real-hyperbolic space-form value
\(h_{\mathrm{cap}}=(n-1)\sqrt\kappa\),

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
process-time gauge \(c=1\). In Lean, the logical separation is explicit:

```text
Addressable β c hcap
  + CapacitySaturated β c hcap
  + hcap_eq_spaceForm hcap n κ
  → normalized state equation
```

Neither saturation nor the space-form identification is part of `Addressable`.
The predicate `hcap_eq_spaceForm` is an algebraic identification of rates, not
a theorem that the host is \(\mathbb H^n_\kappa\).

## Files

```text
ActiveGeometry/
├── Packing.lean         # metric kernel; addressability_limit
├── Capacity.lean        # algebra of the bound (floor, gauge, η)
├── StateEquation.lean   # optional face: saturation + space-form chart
└── Measurability.lean   # growth-class gate identities (instrument)
```

### `Packing.lean`

This file uses Mathlib's canonical `Metric.packingNumber`; it does not
axiomatize a capacity envelope. Formalized results include:

| Declaration | Meaning |
|---|---|
| `card_le_packingCount` | every finite separated subset of a ball is bounded by its exact packing number |
| `subball_fraction_le_packing_fraction` | the fraction of codewords in any smaller sub-ball is bounded by its packing fraction (finite-radius radial concentration) |
| `exists_optimal_blockCode` | an exact finite packing code exists whenever the ball packing number is finite |
| `exists_optimal_blockCode_of_properSpace` | exact finite-block achievability in every proper metric host |
| `hasFinitePacking_of_properSpace` | in any proper metric space the finiteness hypothesis is a theorem |
| `FaithfulRepresentation` | finite separated census in a ball; no retention or asymptotic premise |
| `represented_card_le_packingCount` | represented histories obey the packing bound at every depth |
| `representedRate_le_capacity_eventually` | finite counts induce the normalized rate inequality |
| `addressability_limit_of_hasFinitePacking` | diverging radii and three independent ordinary limits prove `Addressable β c hpack` |
| `addressability_limit` | the same theorem with packing finiteness discharged in a proper host |
| `no_positive_growth_at_zero_capacity` | zero packing capacity excludes positive represented growth |
| `RetainedRepresentation` | faithful representation plus `points_monotone` |
| `represented_card_mono` | retention makes represented counts nondecreasing in depth |

The formal theorem uses ordinary finite limits for represented growth, radial
rate, and packing growth. The full spine's limsup version is a more general
paper theorem. `HasFinitePacking` is a hypothesis of the general theorem and a
proved consequence of `ProperSpace` (ℝⁿ, hyperbolic space, and every complete
Riemannian manifold via Hopf–Rinow), so the intended host class needs no extra
assumption.

The upper bound and `exists_optimal_blockCode` together identify operational
finite-block address capacity exactly with metric packing capacity. This
achievability result is fully metric and host-agnostic, but deliberately does
not assert that optimal codebooks at successive radii are nested, causal, or
preserve a source hierarchy's relational metric. The asymptotic block identity
\(C_{\rm block}(c,\varepsilon)=c h_{\rm pack}\) and stronger constrained
achievability problems remain paper-level statements.

### `Capacity.lean`

Formalized results include:

| Declaration | Meaning |
|---|---|
| `addressability_forces_positive_entropy` | \(\beta>0\), \(c>0\), and \(\beta\le c h_{\rm cap}\) imply \(h_{\rm cap}>0\) |
| `efficiency_le_one` | \(\eta=\beta/(c h_{\rm cap})\le1\) |
| `curvature_at_least_floor` | a space-form capacity inequality gives a curvature lower bound |
| `addressable_spaceForm_floor` | composes `Addressable` with `hcap_eq_spaceForm` |
| `saturated_curvature_eq_floor` | saturation fixes raw curvature once \(c,n\) are fixed |
| `saturated_spaceForm_eq_floor` | composes saturation with a space-form identification |
| `floor_saturates_capacity` | the curvature floor realizes equality |
| `normalized_floor_eq_saturation` | multiplying by \(c^2\) gives the saturation-face value |
| `process_time_gauge` | \(c=1\) recovers the familiar formula |
| `normalized_curvature_scale_invariant` | \(c^2\kappa\) is invariant under radial rescaling |
| `normalizedCurvatureAtSaturation_anti_n` | saturation-face value decreases algebraically with dimension; not a selection of \(n=2\) |

### `StateEquation.lean`

The optional equality-case face. Formalized results include:

| Declaration | Meaning |
|---|---|
| `normalized_state_equation` | saturation plus `hcap_eq_spaceForm` yields the normalized equality |
| `rateMismatchSq_eq_scaled_sqrtMismatch` | the two gap diagnostics are one object up to \((c(n-1))^2\) |
| `rateMismatchSq_zero_iff` | the squared rate gap vanishes exactly on the face |
| `rateMismatch_zero_at_floor` | space-form capacity at the curvature floor equals demand |

These diagnostics are not a Lyapunov theorem and not an evolution law.

### `Measurability.lean`

| Declaration | Meaning |
|---|---|
| `midpoint_exponent_eq` | endpoint-matched log-gap at \(\sqrt{r}\) (Lemma 1.2) |
| `spanInformation_eq_logGap_logMean` | maximum gap is \(\Delta(r,d)\) at the logarithmic mean |
| `logGapDeriv_logMean` | formal derivative vanishes at that mean |
| `spanInformation_pos` | \(\Delta(r,d)>0\) for \(r>1\), \(d>0\) |

Lean does not formalize Hellinger distance, Le Cam's lemma, or the
Poisson-increment model.

## What Lean does not establish

The formalization intentionally does not claim to machine-check:

1. the full limsup generalization of the convergent-rate packing theorem;
2. equivalence of packing and volume entropy under bounded geometry;
3. the space-form classification or the hyperbolic volume formula;
4. Theorem 4.4 (Skenderi / weighted relational capacity of \(\mathbb H^n_\kappa\));
5. Theorem 7.1 (Heintze isotropy / axiom A3);
6. the Buneman/Gromov tree-classification theorems;
7. Sarkar's low-distortion embedding theorem;
8. nested, causal, or relation-preserving achievability;
9. a physical dynamics toward capacity saturation;
10. empirical membership of any biological or linguistic system;
11. alphabet or DNA entropy-rate ceilings (those are substrate constants,
    not the addressability kernel).

These are respectively paper proofs, classical cited results, open modeling
choices, or empirical questions.

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
