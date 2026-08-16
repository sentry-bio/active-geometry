# Paper II outline (v5) — *Evolution as Active Geometry*

Draft for review, not a manuscript. Paper I carries the machinery. This
paper cites it. The v4 outline derived host class and reticulation from
the ladder. v5 is the same derivation, said once.

**The one distinction.** A record is not the thing recorded. It can fail
its original in exactly two ways:

- **Poverty.** The record goes quiet. Splits become unreadable because the
  tape has been overwritten too many times. Recoverable: Test A rises
  through \(L=400{,}000\) with no flattening; decided-quartet accuracy
  extrapolates to an asymptote at or above 1. Nothing is lost in
  principle; you need more tape.
- **Mixing.** The record goes wrong. A second genealogy is written on the
  same tape. Not recoverable: Test C is flat across a \(16\times\) increase
  in \(L\). The transferred block resolves to the donor's topology at
  \(98.7\%\) — the fidelity an untouched block shows for the recipient.
  Zero undecidable quartets at any transfer rate. It does not get
  ambiguous. It gets confidently misdirected.

Those are not two findings. They are one dichotomy, and they organize the
paper:

- **§2** — the two failure modes, on ground truth
- **§3** — a representation that escapes poverty (functional, not sequence)
  still looks room-rich
- **§4** — mixing leaves a scar in real data
- **§5** — failure to distinguish poverty from mixing is why \(\kappa\)
  never stabilized

**Inclusion rule.** A result enters the spine only if it measures this
distinction (quiet vs wrong; recovers with \(L\) vs floors; packing vs
splits) on ground truth or on an independent biological representation.
Simulations license the dichotomy; they do not place a clade. Absolute
curvature, back-solved \(n\), and the state equation do not enter the
headline. Channel saturation (\(p\to 0.75\)) is never reported as
host-utilization saturation (\(\eta\to 1\)).

**Two holds, stated up front.**

1. The dichotomy is verified only in simulation. Test B (misspecification
   alone, no reticulation) is unrun. Until it runs, a non-recovering floor
   on real data cannot be attributed to mixing rather than to a wrong
   substitution model. That is the gate on the biological application.
2. “Genealogy is the expensive thing” is interpretation, not measurement.
   What is measured: one failure recovers and one does not, on synthetic
   data with exact ground truth, plus two independent HGT signatures in
   real data. The framing is a reading of those facts, and is labeled as
   one.

---

## Title (proposed)

*Evolution as Active Geometry: When the Record Fails the Tree*

No state equation, no curvature value, no atlas, no clock. Host class and
reticulation remain the two live tests; they are now named as the two
ways a record fails.

## What Paper I licenses (one paragraph, not a section)

The bound: retained novelty cannot outgrow addressable room. The ladder:
distinguishable endpoints (block) are not preserved genealogy (relational).
The classifier: quartets decide tree-ness; undecidable \(\neq\) confidently
wrong; \(\delta\) is orthogonal to packing. The chart: given A3, real
\(\mathbb H^{d+1}\) inside the Heintze class; \(n=2\) is embeddability.
This paper tests the ladder on biological records. It does not fill the
chart and does not saturate it.

## Abstract — beats

1. A genome is a record of a genealogy, not the genealogy. Paper I's
   ladder says identity and origin are different capacities. This paper
   asks how the record fails the original.
2. On ground truth, it fails in two ways. Poverty: quartets go quiet and
   recover with more sites. Mixing: quartets go confidently wrong and do
   not recover; block distinguishability stays high; a transferred block
   is the donor's tree.
3. A functional representation that is not the sequence tape still looks
   room-rich against a Euclidean null (KEGG; E9 candidate, not E9).
4. Mixing leaves a scar: HGT-enriched functions need extra dimension;
   deployed placement failures are \(7.7\times\) HGT-enriched.
5. Years of unstable \(\kappa\) are what you get if you average quiet
   records with wrong records. No curvature constant is claimed.

## §1 Introduction

- A record is not the thing recorded. Cite Paper I for the ladder, not
  for a curvature of life.
- The record can fail quietly or wrongly. Those demand opposite
  responses — more data versus a different model of the object.
- Roadmap: §2 the dichotomy on ground truth; §3 escaping poverty; §4 the
  mixing scar; §5 why \(\kappa\) never settled; §6 controls; §7
  discussion.

## §2 The two failure modes, on ground truth

Not a biological result. The license.

- **§2.1 Poverty (Test A).** JC69, exact Yule genealogy. Block witness
  count stays near the ceiling. Decided-quartet accuracy rises with \(L\)
  at \(D=2\) and \(D=3\) through \(400{,}000\) sites, no flattening.
  Undecidable fraction shrinks on a \(\ln L\) schedule. Infinite-sites
  at matched event count does not collapse. Mechanism: site reuse, not
  “more mutations.”
- **§2.2 Mixing (Test C).** Non-saturating depth, HGT rate swept. Block
  witness count drifts *up* (\(247\to 254\)). Undecidable fraction is
  identically zero. Decided accuracy floors under \(16\times L\).
  Transferred block vs donor topology \(98.7\%\); vs recipient \(52.3\%\).
- **§2.3 The diagnostic.** \(L\)-recovery: if decided accuracy climbs,
  the record is poor; if it sits, the record is mixed — *provided* Test B
  has ruled out misspecification. Status: reported local simulations;
  reproduce in-repo before citing tables in a submitted draft.

## §3 A representation that escapes poverty still looks room-rich

Host class, derived: if the failure mode in §2.1 is *the tape*, a
non-tape representation should still show exponential room.

- **§3.1** Sequence-native distances at depth are the poor record. Do not
  read host class off them.
- **§3.2** Native functional geometry (KEGG, 1,935 genomes, three
  domains, no geometric model): \(\delta\)-ratio \(0.32\) vs
  dimension-matched Euclidean null. Candidate for E9, not a
  matched-capacity run. Name the gap.
- **§3.3** Synthetic E9-shaped pilot (stress \(0.143\) vs \(0.322\)):
  the test can discriminate. Not biology.
- **§3.4** Not a \(\kappa\), not A3, not saturation. Polar split of an
  imposed embedding is definitional; that trap is why E9 exists.

## §4 Mixing leaves a scar in real data

Classifier, derived: if §2.2 is a second tree on the same tape, real HGT
should show extra angular room and operational failure, without packing
having to fall.

- **§4.1** Prediction: block distinguishability holds; defect does not
  recover with more characters; extra dimension is a readout, not a
  fitted \(n\).
- **§4.2** Admissible now: COG elbows (V/W/B/X, HGT-enriched);
  Atlas placement failures \(7.7\times\) HGT-enriched, zero core-metabolic
  pathway failures.
- **§4.3** recA \(n=3.02\) held out until re-derived without back-solving.
- **§4.4** Next measurement, not this draft: \(L\)-recovery on a real
  deep alignment, after Test B.

## §5 Why \(\kappa\) never stabilized

The archival payoff. Not a curvature section in disguise.

- Years of distance-matrix \(\kappa\) mixed poor records with mixed
  records. Those demand opposite responses. Averaging them is a
  reasonable account of non-convergence.
- Two estimators on the same trees still disagree by \(\sim 2.8\times\)
  (ball-growth vs stress-fit). Named as an open instrument question
  (ladder-rung mismatch). Not offered as anisotropy, not as saturation,
  not as the state equation.
- Comparative leftovers, if kept at all, stay here and stay comparative:
  viral relative curvature tracks depth not mutation rate; cross-alphabet
  ratio as a ratio. No certified absolute \(\kappa\).
- IIb kill lines belong here as the open bet this paper does not take.
  Explicit non-claim: no \(\eta\to 1\).

## §6 Controls

Estimator tests.

- Euclidean null; synthetic \(\kappa\) recovery as meter recovery, not
  biology saturating; destroyed-structure Procrustes; gauge control
  (\(6\times\) swing at fixed \(h_{\mathrm{eff}}\)).
- Minimal encoder: imposed \(\mathbb H^2\), residual \(0.020\), gauge
  \(O(2)\). Stability inside the imposed model, not host-class selection.

A3 / sector-wise packing remains Paper I's empirical remainder,
downstream of E9. One sentence here or in §7, not a competing spine.

## §7 Discussion — the throughline

Rewrite of v4's scope list. Same claims. One argument.

**The reading, labeled as interpretation.** Life keeps two things in one
place: what a thing is (a distinguishable address) and where it came from
(a genealogy). Paper I's ladder says these are different capacities. The
evidence says they degrade independently. Block capacity is robust:
genomes stay distinguishable through saturation and through transfer;
Test C's witness count even rises, because mixing makes sequences more
distinct from their vertical neighbors, not less. Relational capacity is
fragile in two separable ways. A phylogeny is not a record read off the
genome. It is a reconstruction of a physical quantity with a measurable
ceiling, hit either by erosion (fixable) or by structural contradiction
(not). Reticulation is a second tree competing for the same tape, not
noise in a first tree.

**Claimed (I/IIa).** The dichotomy on ground truth; a non-tape
representation still looking room-rich against a Euclidean null
(candidate); two HGT scars in real data, consistent with mixing not
poverty.

**Not claimed.** That Test B is done; that a real-data floor is mixing
rather than misspecification; that “genealogy is expensive” is a
measurement; A3; saturation; a filled atlas; CCS as host-class proof; any
\(\kappa\) constant.

**Falsifiability.** (i) E9 matched-capacity fails on independent real
data → host class fails; (ii) HGT instruments disagree in sign once recA
is re-derived → mixing scar fails; (iii) \(L\)-recovery on a known-HGT
alignment climbs to 1 → the floor prediction fails for that system;
(iv) Test B floors without HGT → the diagnostic does not yet distinguish
mixing from a wrong model.

**Close, not with “biology is active geometry.”** That is the IIb claim
this paper brackets. Close with:

> Life's record keeps identity cheaply and genealogy expensively. When
> the record fails its original, it fails in one of two ways — quietly,
> which more data repairs, or wrongly, which it does not. The distinction
> is measurable, and it is what a geometry of hierarchy is for.

That is Layer I/IIa, it is what the evidence supports, and it leaves the
state equation unearned.

## Sequencing

1. Reproduce A/C in-repo, or cite as reported-not-reproduced.
2. Re-derive recA without back-solving, or drop it.
3. Write §4 from the two admissible HGT instruments.
4. Test B before any real-data floor is called mixing.
5. E9 and the sector test are follow-ups, not prerequisites for the
   dichotomy paper.

## Relative to v4

v4 had the right pieces and too many spines. v5 keeps every claim and
lets the poverty/mixing cut do the organizing. §5 is no longer a
curvature attic; it is the diagnosis of the archive. §7 is the
throughline instead of a not-claimed list (the list remains, inside the
argument).
