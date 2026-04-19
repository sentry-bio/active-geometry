# Formal Proofs in Lean 4

Machine-checked proofs of the core mathematical claims of the geometric
state equation. This Lean stack is the formal backbone of the Hyperbolic
Trilogy:

- **Paper I** — Fenn & Fenn (2026), *Evolution as Active Geometry: The
  Geometric State Equation of the Tree of Life*, bioRxiv
  [10.64898/2026.03.09.710612](https://www.biorxiv.org/content/10.64898/2026.03.09.710612v2) (this repository; DNA/RNA/protein substrate)
- **Paper II** — Fenn & Fenn (2026), *A Geometric State Equation for
  Information-Generating Hierarchies*, Zenodo
  [10.5281/zenodo.19381558](https://doi.org/10.5281/zenodo.19381558)
  (theoretical core; this file ships with it)
- **Paper III** — Fenn & Fenn (in prep), *Convergent Alphabets* ([github.com/sentry-bio/convergent-alphabets](https://github.com/sentry-bio/convergent-alphabets); phoneme substrate)

Part I and Part III of the Lean file are substrate-agnostic (they prove
properties of the state equation for any `h` and `n`). Part II contains
both a DNA-specific block (historical, matches Paper I §3.2) and a
trilogy-general block (`H_raw_of_alphabet α`, `kappa_bounded_by_alphabet_general`)
that covers proteins, phonemes, and any future substrate.

## Layout

```
ActiveGeometry/
└── KappaCurvature.lean   # Full formalization (~370 lines, 0 sorries)
```

The file is organized in three parts:

- **Part I — Core state equation.** The curvature formula
  `κ(h, n) = (h · ln 2 / (n − 1))²`, its specialization at `n = 2`,
  positivity, uniqueness, monotonicity in `h` and `n`, maximization at
  `n = 2`, and the growth-rate matching identity `r · √κ = h · r · ln 2`.

- **Part II — Entropy rate decomposition.** The four-letter channel
  capacity `H_raw := log₂ 4 = 2` (matching paper §3.2), the lower anchor
  `h_three_letter := log₂ 3 ≈ 1.585` with a proof that `1.5 < log₂ 3 < 1.7`,
  a multiplicative decomposition `h_effective = H_raw · φ · ψ · ω` for
  transition/context/selection biases in `(0, 1]`, and the resulting
  alphabet-capacity bound `κ(h_eff, 2) ≤ 4·(ln 2)² ≈ 1.921`.

- **Part III — Tree dimensionality and Lyapunov stability.** The
  four-point condition definition of a metric tree, the rate-distortion
  residual `ε = h · ln 2 − (n − 1)√κ`, the potential `U = ε²` with its
  unique zero at `κ_critical`, and the Lyapunov function
  `V(κ, κ*) = (√κ − √κ*)²` establishing global stability of the state
  equation solution.

## Building

```bash
cd theory/lean
lake build
```

Requires Lean 4 via `elan` and a Mathlib 4 cache.

## Theorem Inventory

| Theorem | Statement |
|---|---|
| `kappa_n2` | `κ(h, 2) = (h · ln 2)²` |
| `kappa_pos` | `h > 0 ⟹ κ(h, 2) > 0` |
| `kappa_unique` | unique positive root of the state equation |
| `kappa_mono_h` | `h₁ < h₂ ⟹ κ(h₁, 2) < κ(h₂, 2)` |
| `kappa_mono_n` | `n₁ < n₂ ⟹ κ(h, n₂) < κ(h, n₁)` |
| `kappa_max_at_n2` | `n > 2 ⟹ κ(h, n) < κ(h, 2)` |
| `kappa_scaling` | `κ(c·h, 2) = c² · κ(h, 2)` |
| `growth_rate_match` | `r · √(κ(h,2)) = h · r · ln 2` |
| `H_raw_eq_two` | `H_raw = 2` (four-letter channel capacity) |
| `h_three_letter_bounds` | `1.5 < log₂ 3 < 1.7` |
| `h_three_letter_lt_H_raw` | `log₂ 3 < H_raw = 2` |
| `entropy_rate_decomposition_bounds` | `0 < h_eff ≤ H_raw` |
| `kappa_bounded_by_raw` | `κ(h_eff, 2) ≤ κ(H_raw, 2)` |
| `kappa_bounded_by_alphabet` | `κ(h_eff, 2) ≤ 4·(ln 2)²` (DNA, α = 4) |
| `H_raw_of_alphabet` | `α ↦ log₂ α` (trilogy-general capacity) |
| `H_raw_of_alphabet_pos` | `α > 1 ⟹ H_raw_of_alphabet α > 0` |
| `H_raw_of_alphabet_mono` | strictly monotone in α |
| `H_raw_eq_alphabet_four` | `H_raw = H_raw_of_alphabet 4` |
| `kappa_bounded_by_alphabet_general` | `κ(h_eff, 2) ≤ (log α)²` for any α > 1 |
| `potential_nonneg` | `U(h, κ, n) ≥ 0` |
| `potential_zero_iff` | `U = 0 ⟺ κ = κ_critical` |
| `potential_gradient_zero_at_critical` | `ε(h, κ_critical, n) = 0` |
| `lyapunov_nonneg` | `V(κ, κ*) ≥ 0` |
| `lyapunov_zero_iff` | `V = 0 ⟺ κ = κ*` |
| `embedding_dimension_optimal` | `n > 2 ⟹ κ(h, n) < κ(h, 2)` |
| `dimension_two_maximizes_curvature` | `∀ n > 2, κ(h, n) < κ(h, 2)` |

All theorems are proved. No `sorry` and no extra axioms beyond
`propext`, `Classical.choice`, and `Quot.sound` (standard Lean/Mathlib).

## What Is Not in Lean

Two claims live outside the Lean file because they are physical, not
mathematical:

1. **`n = 2` as the embedding dimension of life.** The state equation
   admits any `n ≥ 1`. That evolution sets `n = 2.00 ± 0.05` across all
   tested systems is an empirical finding (paper §2, Table 1), not a
   theorem.

2. **`h ≈ 1.61` bits as the entropy rate of DNA.** The decomposition
   `h_effective = H_raw · φ · ψ · ω` is formalized with biases in
   `(0, 1]`, but the specific values of the three biases (transition bias,
   CpG context, purifying selection) are empirical inputs from molecular
   biology (paper §3.2).

Feed an independently measured `(h, n)` into the state equation and the
Lean proofs guarantee the resulting `κ` is unique, positive, and
Lyapunov-stable. The measurement is biology's job; the machinery above
makes sure the math is honest once the measurement is made.

## References

1. Fenn, R. & Fenn, A. (2026). *Evolution as Active Geometry: The
   Geometric State Equation of the Tree of Life.* bioRxiv 2026.03.09.710612.
2. Manning, A. (1979). Topological entropy for geodesic flows.
   *Ann. Math.* **110**, 567–573.
3. Sarkar, R. (2012). Low-distortion Delaunay embedding of trees in the
   hyperbolic plane. *Graph Drawing 2011*, LNCS **7034**, 355–366.
4. Gromov, M. (1987). Hyperbolic groups. *Essays in Group Theory*.
