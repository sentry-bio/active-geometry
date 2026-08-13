/-
  Active Geometry: Main Entry Point
  ==================================

  The coordinate-free mathematical claim is the addressability bound

    β ≤ c · h_cap,

  where β is retained-information growth in nats per generative step,
  c converts generative steps to radial distance, and h_cap is host packing
  entropy (or volume entropy under additional hypotheses). Capacity saturation
  and isotropic realization are separate predicates. For β > 0, c > 0, n > 1,
  and κ ≥ 0, an isotropic hyperbolic host implies

    κ ≥ (β / (c · (n - 1)))².

  Capacity saturation gives equality. For β = h · ln 2 and normalized
  sectional-curvature magnitude κ̄ = c²κ:

    κ̄ = (h · ln 2 / (n - 1))².

  The familiar raw-curvature-magnitude formula uses the process-time gauge
  c = 1.
  See theory/MATHEMATICAL_SPINE.md for hypotheses, proofs, scope, and
  the independent role of the four-point tree condition.
-/

import ActiveGeometry.Addressability
import ActiveGeometry.Packing
import ActiveGeometry.KappaCurvature
