# Paper II outline (v4) — *Evolution as Active Geometry*

Draft for review, not a manuscript. Derived from the v3 outline, the
program map ([`PROGRAM.md`](PROGRAM.md), [`THROUGHLINE.md`](THROUGHLINE.md)),
and the synthetic A/C simulations (reported locally; code not yet in this
repository). Paper I carries the machinery. This paper cites it and tests
what it licenses on life.

**The one-line change from v3.** v3 already leads with host class and
reticulation rather than curvature. v4 makes those two questions *derived*:
they are the two axes of Paper I's empirical figure (growth class × tree
defect) applied to a representation that is not the inferred tree. The
simulations do not place a clade. They license the instruments and stop a
misreading: multiple-hit saturation is resolvability poverty; horizontal
transfer is a structural floor. Paper II's biology is organized around that
distinction, not around \(\kappa\).

**Inclusion rule.** A result enters the spine only if it is dimensionless
(a ratio, a correlation, a comparative statement, a null-rejection) **and**
it answers host class or tree defect on an independent representation. A
simulation licenses a prediction and an instrument, not a biological
finding. Absolute curvature, back-solved \(n\), and state-equation
explanations stay in §5, flagged. Channel saturation (\(p\to 0.75\)) is
never reported as host-utilization saturation (\(\eta\to 1\)).

---

## Title (proposed)

*Evolution as Active Geometry: Host Class and Reticulation in the Tree of Life*

No state equation, no curvature value, no atlas, no clock.

## What Paper I licenses (do not re-prove)

Cite once, in this order, because it is the derivation of the paper:

1. **Limit.** \(\beta\le c\,h_{\mathrm{pack}}\). Polynomial hosts cannot
   retain exponential novelty at finite rate (Corollary 4.3).
2. **Ladder.** Block capacity (distinguishable endpoints) is not relational
   capacity (preserved genealogy). Slack splits into tax \(\Gamma\) and
   utilization. In a weighted hyperbolic host, \(\Gamma=0\) as a supremum
   (Theorem 4.4). On a sequence metric that need not hold.
3. **Classifier.** Quartets decide tree-ness, not curvature (Theorem 6.1).
   Undecidable \(\neq\) confidently wrong. \(\delta\) is orthogonal to packing.
4. **Chart.** Given isotropy, space-form plus packing select
   \(\mathbb H_\kappa^n\). Theorem 7.1: inside the Heintze class, axiom A3
   (full \(O(d)\) directional symmetry) forces real hyperbolic space. \(n=2\)
   is embeddability, not a fit. Polar radius/angle are encoder choices.

This paper tests (1)–(3) on biological representations. It does not claim
(4)'s A3, does not fill the chart, and does not saturate it.

The honest figure is the phase diagram already in
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md):

| | tree-like \(\delta\approx 0\) | reticulate \(\delta>0\) |
|---|---|---|
| **exponential room** | hyperbolic-tree candidate | mixed / HGT-rich host |
| **polynomial room** | path-like | networked and capacity-limited |

Only the upper-left cell invites the chart. Sequence distances at depth can
land in the upper-right cell for two different reasons. The simulations
tell them apart. That is the biological application.

## Abstract — beats

1. Paper I forces a room-rich host for any process that retains a
   distinguishable genealogy. This paper asks the two questions that figure
   actually contains: does a *representation* of life (not the inferred tree
   read as its own metric) supply exponential room, and does reticulation
   register as a structural defect rather than as noise.
2. Native functional geometry (KEGG orthology, no geometric model) is
   negatively curved against a dimension-matched Euclidean null — a
   real-data *candidate* for the host-class test, not a run of the
   matched-capacity protocol.
3. Reticulation has a synthetic mechanism and a biological signature.
   Synthetically: JC multiple hits make quartets *undecidable* and recover
   with more sites; HGT makes quartets *confident and partly wrong* and does
   not recover with sites, while block distinguishability stays high. In
   data: HGT-enriched functions need extra predictive dimension; deployed
   placement failures are HGT-enriched. Excess dimension is a readout of
   extra angular room, not a back-solved \(n\).
4. No curvature constant is claimed. Two estimators of packing rate on the
   same trees disagree by \(\sim 2.8\times\); reported as an open instrument
   question (ladder-rung mismatch), not as anisotropy and not as
   saturation.
5. Whether the host is isotropic is axiom A3, empirical, unrun. The test is
   sector-wise packing after host class is known.

## §1 Introduction

- One paragraph: remembering while creating costs room. Cite Paper I. This
  paper does not prove the bound.
- Second paragraph: a finite-alphabet *record* of a tree is not the tree.
  Endpoints can stay distinct while splits go quiet (channel poverty) or
  while a second tree is written on the same tape (transfer). That is the
  ladder, and it is why host class and reticulation are different sections.
- Third: no curvature of life is claimed. Roadmap: §2 derivation from
  simulation (instruments, not clades); §3 host class; §4 reticulation;
  §5 curvature bounded; §6 A3; §7 controls; §8 scope.

## §2 What the simulations license (not a biological result)

New relative to v3. Short. Grounds the application without smuggling
synthetic numbers into the spine.

- **§2.1 Two failures of relational capacity, one packing.** On Yule trees
  with exact genealogy, JC69 keeps almost all tips distinguishable while
  quartet resolvability collapses with depth. Infinite-sites, matched on
  event count, does not collapse. So the binding constraint on a
  four-letter tape is relational, not block. Mutual exhaustion did not
  appear at realistic \(L\) vs \(N\).
- **§2.2 The \(L\)-recovery test.** Decided-quartet accuracy under multiple
  hits *rises* toward 1 as \(L\) grows to \(4\times 10^5\) (Test A):
  information dilution, an instrument limit, not a host-class failure.
  Under HGT at non-saturating depth (Test C), decided accuracy *floors*
  under a \(16\times\) increase in \(L\); every quartet is computable; a
  growing fraction is wrong; block witness count drifts up. The transferred
  block agrees with the donor topology, not the recipient. That is
  structural \(\delta\), orthogonal to packing — Paper I's E7 prediction,
  instantiated.
- **§2.3 What this applies to biology.** Do not read host class off
  saturated marker-gene distances (they can look non-tree because they are
  poor). Do not call channel saturation \(\eta\to 1\). Do use \(L\)-recovery,
  or an independent representation, to tell poverty from mixing. The
  biological sections below are those two tests, as far as present data
  allow. Status: reported local simulations; reproduce in-repo before citing
  numerical tables in a submitted draft.

## §3 Is the native representation room-constrained? (Layer IIa)

v3's §2, now after the derivation.

- **§3.1 Three geometries.** Imposed, learned, native. Only native distances
  test host class. Sequence-native at large \(D\) is the left branch of the
  throughline and may be poor (§2.2); functional-native (KEGG) is the
  better E9-shaped input because it is not the inferred tree and not JC
  distance.
- **§3.2 Native functional geometry.** KEGG orthology, 1,935 genomes, three
  domains, no geometric model. \(\delta\)-ratio \(0.32\) vs
  dimension-matched Euclidean null; distances monotone across ranks. State
  as a first real-data *candidate* for E9, and name the gap: not
  matched-capacity, not the certified meter protocol.
- **§3.3 Synthetic E9-shaped pilot.** Hyperbolic vs Euclidean, matched
  dimension, same optimizer: stress \(0.143\) vs \(0.322\). Methods
  validation that the test can discriminate, not biology.
- **§3.4 What this does not establish.** Not a \(\kappa\), not A3, not E9
  proper, not saturation. Polar split of an imposed embedding is
  definitional (Euclidean control separated as well or better in the
  Line A polar attempt); that trap is why E9 exists.

## §4 Reticulation as structural defect, not noise (classifier)

v3's §3, re-derived from §2.2 rather than from back-solved \(n\).

- **§4.1 The prediction, from Paper I plus the simulation.** A second tree
  on the same tape should: keep block distinguishability; drive quartet
  defect that does **not** recover with more characters; optionally demand
  extra angular dimension. Packing rate need not fall (E7:
  \(\delta\perp h_{\mathrm{pack}}\)). Excess PCA dimension is a *readout*
  of extra angular room, not the definition of reticulation and not a
  fitted \(n\).
- **§4.2 Primary biological evidence (admissible now).**
  - COG predictive-dimensionality elbow: V/W/B/X need full dimension;
    HGT-enriched; mean elbow \(k=57.7\) vs \(129\) at the tail. No
    \(\kappa\), no state equation.
  - Atlas placement failures \(7.7\times\) HGT-enriched; zero core-metabolic
    pathway failures. Operational, independent.
- **§4.3 recA — hold out until re-derived.** \(n_{\mathrm{backsolved}}=3.02\)
  is inadmissible (n=2 seam). Blocking task: stress vs embedding dimension
  with no \(\kappa\) or \(h_{\mathrm{eff}}\). Until then, two-of-three
  instruments.
- **§4.4 Derived next measurement (not required for v3 submission).** The
  \(L\)-recovery test on a real deep alignment with an independent reference
  topology, plus a misspecification control (Test B, unrun). That would be
  E7's synthetic instrument applied in vivo. Do not pad the paper with a
  rushed version. Name it as the clean follow-up.
- **§4.5 Synthesis.** Strongest finding because it needs no curvature
  number: mixing registers as extra angular room and as operational
  failure, consistent with a relational floor, not with channel poverty.

## §5 Curvature, honestly bounded (Layer IIb, open)

v3's §4. Keep the blanket rule. Add the two-saturations sentence.

- **§5.1** No certified \(\kappa\). Every reported value carries isotropy
  premise and estimator identity.
- **§5.2** Viral relative curvature tracks depth (\(\rho=0.84\)), not
  mutation rate (\(\rho=0.12\)). The \(r=0.996\) entropy correlation is
  secondary and is **not** the state equation (would need A3 and
  \(\eta=1\)).
- **§5.3** Cross-alphabet ratio \(3.1\times\) as comparative, not as two
  certified absolutes.
- **§5.4** Estimator gap \(2.84\times\) (ball-growth vs stress-fit) as open
  instrument question. Named hypothesis: ladder-rung mismatch (packing
  occupancy \(\neq\) weighted-clock relational code). Do **not** offer this
  as A3 failure; that would confound instrument with host. Both remain
  untested; say so.
- **§5.5** Explicit non-claim: no saturation, no \(\eta\to 1\). Independent
  IIb kill lines (domain \(r=0.35\), protein \(-0.11\), viral
  \(0.06\)–\(0.19\)) stay in this section as the open bet, not a buried
  negative. Under a hyperbolic-host assumption those gaps are utilization;
  on a sequence metric they may be tax. This paper does not average them
  into a \(\kappa\).

## §6 The open axiom (A3)

v3's §5, tightened.

- **§6.1** Cite Theorem 7.1. Do not re-derive Heintze. Geometry alone does
  not select the isotropic member; A3 does, as a theorem, inside the
  Heintze class.
- **§6.2** Translate A3 into biology without a kitchen sink. The simulation
  that actually violated directional exchangeability is **structured HGT**
  (a preferred donor). Substitution-spectrum biases (ti/tv, CpG) are
  *channel* properties and may look like Test A, not like a Heintze pinch.
  Do not lump them with A3 until a sector test says so.
- **§6.3** Decisive test, unrun, downstream of E9: sector-wise
  \(h_{\mathrm{pack}}\) vs shuffle-angle null. Companion note or v4. Not a
  rushed figure here.

## §7 Controls

Estimator tests, not curvature claims.

- Euclidean null: \(\kappa=0\) on polynomial random trees, 5/5.
- Synthetic recovery of known \(\kappa_{\mathrm{true}}=(\ln b)^2\) to
  \(1.08\%\) — report as meter recovery, not as biology saturating.
- Destroyed structure: Procrustes \(0.94\to <0.3\).
- Minimal encoder: imposed \(\mathbb H^2\), \(\kappa\) frozen, residual
  \(0.020\), leftover gauge \(O(2)\). Reproducibility **within the imposed
  model**, not host-class selection (that is E9 / §3).
- Gauge control: branching factor fixed, \(\mu\) swept, measured rate
  swings \(6\times\), invariant \(\mathrm{rate}\times\delta\times\ln 2
  \approx \ln b\). Validates `process_time_gauge`. Any rate in this paper
  without a named clock is not comparable across conditions.

## §8 Discussion and scope

- **Claimed.** On independent representations: native functional geometry
  prefers a negatively curved host against a Euclidean null (candidate, not
  E9); reticulation registers as extra angular room and as HGT-enriched
  operational failure, consistent with a relational floor.
- **Licensed but not claimed as biology.** The A/C mechanism: poverty vs
  mixing, diagnosed by \(L\)-recovery.
- **Not claimed.** A curvature constant; A3; saturation; a filled atlas; that
  sequence space at depth is the host; that CCS Procrustes selects the host.
- **Falsifiability.** (i) matched-capacity E9 fails on real independent
  data → host class fails; (ii) HGT instruments disagree in sign once recA
  is re-derived without back-solving → reticulation prediction fails;
  (iii) sector test finds no angular dependence under known directional
  bias → A3-as-stated needs revision; (iv) \(L\)-recovery on a known-HGT
  alignment climbs to 1 → the floor prediction fails for that system.

## Sequencing

1. Reproduce A/C in-repo, or cite them only as “reported, not reproduced
   here,” before any submitted numerical table from them.
2. Re-derive recA without back-solving, or drop it.
3. Write §4 and §5 from existing admissible material.
4. Do not wait on E9 or the sector test to submit a v3-shaped paper; do
   wait on (1)–(2) for honesty. Matched-capacity E9 on GTDB/KEGG remains
   the highest-leverage follow-up, not a prerequisite for the two
   comparative claims.
5. Reconcile viral and protein tables to one provenance (still blocking
   for §5).

## Removed from v3 / added here

| Change | Reason |
|---|---|
| New §2 on simulations | Application to biology is derived from the ladder, not hung on KEGG/HGT after the fact. |
| Reticulation prediction restated as \(\delta\)-floor \(\perp\) packing | Matches E7 and Test C; excess dimension is a readout. |
| Channel poverty vs mixing named | Prevents saturated 16S from being reported as reticulation or as \(\eta\). |
| A3 candidates restricted; HGT primary | Test C is the mechanism; ti/tv is a different rung. |
| Estimator gap *not* offered as anisotropy | Would confound L2a.6 with I.3. |
| Encoder control demoted to imposed-model stability | Throughline / CCS seam. |
| IIb kill lines placed in §5 | Honest open bet, not omitted. |
