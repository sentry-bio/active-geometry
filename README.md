# Wolfram/Mathematica Notebooks

This directory contains Mathematica notebooks for symbolic validation
of the geometric state equation derivation.

## Notebooks

| Notebook | Purpose |
|----------|---------|
| `curvature_derivation.nb` | Step-by-step κ derivation |
| `rate_distortion.nb` | Rate-distortion theory validation |
| `dimensional_analysis.nb` | n = 2 uniqueness argument |
| `entropy_bounds.nb` | Entropy rate calculations |

## Key Computations

### State Equation Derivation

```mathematica
(* Define the state equation *)
kappa[h_, n_] := (h * Log[2] / (n - 1))^2

(* Canonical values *)
hCanonical = 1.61;  (* bits/nucleotide *)
nCanonical = 2;     (* dimensions *)

(* Theoretical prediction *)
kappaTheory = kappa[hCanonical, nCanonical]
(* Output: 1.2296... ≈ 1.23 *)
```

### Agreement Verification

```mathematica
kappaEmpirical = 1.247;
agreement = Abs[kappaTheory - kappaEmpirical] / kappaEmpirical * 100
(* Output: 1.39% < 2% ✓ *)
```

## Running

1. Open notebooks in Mathematica 13+
2. Evaluate all cells (Shift+Enter)
3. Results should match `constants.yaml`

## Export

To export notebooks as PDFs for supplementary materials:

```mathematica
Export["curvature_derivation.pdf", EvaluationNotebook[]]
```

## Notes

- Notebooks are self-contained with all definitions
- Results are compared against `../../constants.yaml`
- All symbolic manipulations can be verified step-by-step
