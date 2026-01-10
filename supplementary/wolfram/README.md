# Wolfram Language Verification (Supplementary)

These Wolfram Language (.wl) files provide Computer Algebra System (CAS) verification
of the mathematical derivations. They are **supplementary** material—the formal proofs
are in Lean 4 (`theory/lean/`) and the numerical validation is in Python (`validation/notebooks/`).

## Files

| File | Purpose |
|------|---------|
| `SI2_Wolfram_Skeleton.wl` | Master CAS verification (consolidates A-D) |
| `NotebookA_SelfConsistency.wl` | Self-consistency and uniqueness proofs |
| `NotebookB_SensitivityRobustness.wl` | Sensitivity analysis, error propagation |
| `NotebookC_CrossDomainPrediction.wl` | Cross-domain kappa predictions |
| `NotebookD_NullSimulations.wl` | Null simulation functions |

## Running

Requires Wolfram Mathematica or Wolfram Engine (commercial license).

```mathematica
(* In Mathematica *)
<< "SI2_Wolfram_Skeleton.wl"
```

## Note

The Python notebooks in `validation/notebooks/` replicate all numerical results
using SymPy and mpmath, requiring no commercial license. The Lean proofs in
`theory/lean/` provide formal machine-checked verification.
