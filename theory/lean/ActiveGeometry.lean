/-
  Active Geometry: public entry point
  ===================================

  Cite `Packing.addressability_limit` (re-exported below as
  `addressability_limit`). That is the convergent-rate packing theorem in
  every proper metric host: a faithful finite-resolution representation with
  convergent represented growth `β`, radial rate `c`, and host packing growth
  `h_pack` satisfies

      β ≤ c · h_pack.

  What this library machine-checks:

  * the packing converse, as ordinary finite limits (`Tendsto`);
  * the finite-block identity `A_block = packing number` in proper spaces;
  * the real algebra of the bound (efficiency, curvature floor, gauge);
  * the equality-case face, *given* saturation and a space-form identification
    of `h_cap`.

  What this library does not machine-check:

  * the paper's limsup generalization of the packing theorem;
  * Theorem 4.4 (Skenderi / weighted relational capacity of `ℍⁿ_κ`);
  * Theorem 7.1 (Heintze isotropy / axiom A3);
  * packing entropy = volume entropy under bounded geometry;
  * nested, causal, or relation-preserving achievability;
  * a dynamics toward saturation, or occupancy of any biological host.

  File cut:

  * `Packing`          — metric kernel; `addressability_limit`
  * `Capacity`         — algebra of the bound (floor, gauge, `η`)
  * `StateEquation`    — optional face: saturation + space-form chart
  * `Measurability`    — growth-class instrument, not a theory layer

  There is no DNA, alphabet, or biosphere-curvature namespace in this library.
  See `theory/ADDRESSABILITY_KERNEL.md` for the paper boundary.
-/

import ActiveGeometry.Capacity
import ActiveGeometry.Packing
import ActiveGeometry.StateEquation
import ActiveGeometry.Measurability

namespace ActiveGeometry

export Packing (addressability_limit)

end ActiveGeometry
