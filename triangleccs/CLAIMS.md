# Claim ledger

Status vocabulary (compatible with Active Geometry theory/CLAIMS.md):

- **THEOREM** — proved elsewhere; cited, not re-derived here.
- **CONVENTION** — gauge / coordinate choice, true by stipulation.
- **INSTRUMENT** — status of a measurement tool.
- **EMPIRICAL** — decided only by experiment.
- **CANDIDATE** — provisional until the freeze-gate passes.
- **ADVISORY** — carried and labelled; not certified.
- **CIRCULAR** — semi-circular; flagged, not leaned on.
- **OVERLAY** — optional Layer IIb / uniqueness material; not form.

| # | Claim | Status | Where |
|---|---|---|---|
| C1 | Retained novelty obeys β ≤ c · h_pack | THEOREM (cite AG Packing.lean) | `docs/WHY.md` |
| C2 | Block address capacity = packing (finite radius) | THEOREM (cite AG) | `triangleccs/packing/bound.py` |
| C3 | Weighted relational capacity of H_κ^n equals block capacity | THEOREM (cite AG Thm 4.4; existence, not uniqueness) | `triangleccs/chart/poincare.py` |
| C4 | Four-point condition ⇔ additive tree metric | THEOREM (cite AG) | `triangleccs/classifier/quartets.py` |
| C5 | κ on Form is a frozen gauge | CONVENTION | `triangleccs/datum/form.py` |
| C6 | dim = 2 is inhabit H² | CONVENTION | `triangleccs/datum/form.py` |
| C7 | Tokenizer is a frame parameter | CONVENTION | `triangleccs/datum/form.py` |
| C8 | θ is the certified-candidate axis | CANDIDATE until freeze-gate | `triangleccs/address.py` |
| C9 | r is advisory | ADVISORY | `triangleccs/address.py`, `triangleccs/chart/polar.py` |
| C10 | Balloon: sequence metric can keep endpoints while losing quartets | INSTRUMENT (Yule/JC69 on-tree) | `triangleccs/tape/balloon.py` |
| C11 | State equation / η → 1 | OVERLAY | overlays/ (empty in v1) |
| C12 | Curvature genericity | OVERLAY | overlays/ (empty in v1) |
| C13 | No public (h, kappa) → n | INSTRUMENT (enforced) | `tests/test_firewall.py` |
| C14 | Genomic LM is not the MDL of diversity; sextant is the chart encoder | CONVENTION (this datum) | `docs/ENCODER.md`, `triangleccs/sextant/place.py` |
| C15 | Chart packing uses the Poincaré metric, not Euclidean coordinates | INSTRUMENT | `triangleccs/metric.py`, `triangleccs/packing/bound.py` |
