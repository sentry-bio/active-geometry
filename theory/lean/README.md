# Formal Proofs in Lean 4

This directory contains formally verified proofs of the theoretical claims
in "A Geometric State Equation for Evolutionary Dynamics".

## Structure

```
ActiveGeometry/
├── Basic.lean       # Foundational definitions
├── Curvature.lean   # Main curvature theorems
└── StateEquation.lean # State equation derivation
```

## Building

```bash
# Install Lean 4 via elan
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh

# Build proofs
lake build
```

## Key Theorems

### Theorem I: State Equation
```lean
theorem state_equation_positive (h : EntropyRate) (n : Dimensionality) :
    0 < stateEquation h n
```

The curvature κ = (h ln 2 / (n-1))² is uniquely determined and positive.

### Theorem II: Dimensionality
```lean
axiom four_point_implies_n2 : ...
```

Phylogenetic trees (4-point condition) embed isometrically in H².

### Theorem III: Agreement
```lean
theorem agreement_bound :
    |kappa_theory - kappa_empirical| / kappa_empirical < 0.02
```

Theory and measurement agree within 2%.

## Dependencies

- Lean 4
- Mathlib4

## Verification Status

| Theorem | Status |
|---------|--------|
| State equation positivity | ✅ Proved |
| Dimensional uniqueness | 📝 Axiomatized |
| Agreement bound | 🔢 Numerical |

## References

1. Fenn & Fenn (2025), "A Geometric State Equation for Evolutionary Dynamics"
2. Manning, A. (1979), "Topological entropy for geodesic flows"
