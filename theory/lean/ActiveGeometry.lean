/-
  Active Geometry: Main Entry Point
  ==================================

  This module contains machine-checked proofs of the geometric
  state equation for evolutionary dynamics:

    κ = (h ln 2 / (n-1))²

  where:
    κ = hyperbolic curvature
    h = entropy rate (bits/symbol)
    n = embedding dimension

  For n = 2: κ = (h ln 2)² ≈ 1.247 when h ≈ 1.61 bits/symbol
-/

import ActiveGeometry.KappaCurvature
