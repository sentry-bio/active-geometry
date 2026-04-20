/-
  Active Geometry: Main Entry Point
  ==================================

  Machine-checked proofs of the geometric state equation for the
  tree of life (Fenn & Fenn 2026, bioRxiv 10.64898/2026.03.09.710612):

    κ = (h · ln 2 / (n - 1))²

  where:
    κ = sectional curvature of the embedding hyperbolic manifold
    h = Shannon entropy rate of the generating code (bits/symbol)
    n = embedding dimension

  At n = 2 (the empirical invariant across DNA, RNA, and protein
  alphabets per paper §2), this reduces to κ = (h · ln 2)². For DNA
  with h ∈ [1.58, 1.65] bits/nt, the predicted κ ∈ [1.20, 1.31]
  (paper §3.4). Encoder-free post-hoc telescope sweeps (§4.2) land
  at κ ≈ 1.28–1.34, within the predicted interval.
-/

import ActiveGeometry.KappaCurvature
