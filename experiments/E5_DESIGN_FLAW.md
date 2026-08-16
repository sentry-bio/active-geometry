# E5 design flaw — report, do not push past

**Status.** Instrument unverified. Do not report \(\eta\) from the current
E5 designs. This is not a tuning problem.

**Question E5 was supposed to ask.** Does description-length pressure, in a
non-biological host, drive utilization toward 1? That requires a
representation that actually has exponential packing capacity, so the
pressure term has something to act on.

**What failed.** The learned embedding is not an exponential-capacity
chart of the tree. At \(\beta_{\mathrm{rate}}=0\) (no pressure) the
growth-class gate reads **polynomial**, not exponential. A converged
hyperbolic tree embedding would already be exponential before any
pressure term is applied. This one is not.

**Diagnostics (pressure-trained embedding, not guessed).**

| Check | Observed | What a real \(\mathbb H^2\) tree chart would show |
|---|---|---|
| Spearman of embedded vs true tree distances | \(0.635\) | well above that (near isometric on hops) |
| Radius vs generation depth | \(0.83\) | high, and not sufficient by itself |
| Distance range | \(337\times\) | span is not the issue |
| Ball-growth from the origin | not exponential | exponential (Sarkar / E1–E2) |

Radius can track depth while the **shape** of occupancy fails. Exponential
room in \(\mathbb H^2\) is angular fan-out: children subdivide the parent's
sector. Relative distance ordering plus a radius penalty never enforces
that. Equal-angle shells (this repo's first runner) and sparse triplets
(the pressure-trained run) are both weak for the same reason.

**Why \(\eta\) from this design is not a measurement.** Utilization is
\(\beta/(c\,h_{\mathrm{pack}})\). If \(h_{\mathrm{pack}}\) is not the packing
rate of an exponential-capacity representation of the generator, the ratio
does not answer E5. There is no reliable exponential-capacity
representation here for the pressure term to act on.

**This repo's first runner** (`experiments/e5_trained_hierarchy.py`,
`experiments/e5_small_scale.json`) used pairwise hop-MSE and equal-angle
shells, not Sarkar sector subdivision. It reported a bake-in kill
(shuffled targets still exponential; \(\lambda\) sweep flat). That kill
was the right *refusal* of \(\eta\). The mechanism above is the right
*diagnosis*: occupancy was never shown to be recovered hyperbolic tree
geometry. Do not read the uncertified \(\eta\approx 0.93\) band as
near-saturation.

**What a redesign must include (either, probably both).**

1. **Sarkar-style initialization / layout** — children subdivide the
   parent's angular sector, as in the E1/E2/E9 hyperbolic tree embedder.
   A fixed cone half-angle is not admissible (E2).
2. **Full pairwise log-distance regression** onto the generator metric,
   not sparse triplets and not ordering-plus-radius alone.

Until an embedding of the *same* generator passes (i) Spearman vs tree in
the isometric regime, (ii) exponential ball-growth from the origin at
\(\beta_{\mathrm{rate}}=0\), and (iii) polynomial growth on a grid control,
E5 does not run. It has an unverified instrument.

**Not claimed.** Genomic saturation. A state equation. That the limit
theory is damaged. Paper II is untouched.
