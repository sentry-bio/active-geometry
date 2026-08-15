# Active Geometry — the program, whole

This is the map of the whole artefact: the layered theory, the two manuscripts
it supports, the seams that join them, and the machinery that makes the whole
self-verifying. It is deliberately short; each part points to the document that
carries the detail.

## The one idea

Remembering while creating costs room. A process that generates and *retains*
distinguishable histories can never gain them faster than its representation
space gains room to keep them apart, and "room" is a counting property of
geometry. Everything else is this inequality, its exact achievable form, the
geometry that realizes it, and the question of whether nature fills it.

## Four layers

The program is stratified so that a failure in one layer cannot contaminate
another. The definitions live in [`CLAIMS.md`](CLAIMS.md); the mathematics in
[`ADDRESSABILITY_KERNEL.md`](ADDRESSABILITY_KERNEL.md) and
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md); the relational coding theorem
in [`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md); the
finite-sample layer in [`MEASURABILITY.md`](MEASURABILITY.md); the tests in
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md).

- **Layer 0 — finite-sample measurability.** What a finite pointed sample
  can decide. The growth-class identities, the Le Cam bound that converts
  radial span into a yes/no, and the instrument that refuses unmeasurable
  matrices. No curvature, no clock, no saturation. This is the only
  mathematical layer whose hypotheses *are* the measurement.

- **Layer I — universal capacity theory (curvature-free).** The bound
  \(\beta\le c\,h_{\mathrm{pack}}\); the exact block identity
  \(C_{\mathrm{block}}=c\,h_{\mathrm{pack}}\); the constrained-capacity ladder;
  the slack decomposition. Proven and machine-checked. Mentions no curvature.

- **Layer IIa — host class (realization).** *Which* geometry hosts the data:
  hyperbolic, at what dimension, and is that *forced*. The space-form argument,
  the curvature-genericity conjecture, \(n=2\) as an embeddability floor. The
  **better-supported** biological claim — but its premise, isotropy, is
  asserted, not measured. Decisive test: E9.

- **Layer IIb — saturation.** *Whether* a process fills its budget, giving the
  state-equation equality. The **harder, less-supported** claim; every
  independent test to date sits below its kill line.

The four-point (Buneman ≡ Gromov) classifier is an independent bridge to IIa:
it decides tree-ness and measures reticulation, but calibrates nothing.

## Two manuscripts

The program publishes as two papers with a firewall between them, so the theory
stands whether or not the biology saturates.

### Paper I — *The Addressability Limit: A Packing Bound for Information-Generating Hierarchies*

- **Domain.** Metric geometry, information theory. No biology in the load-bearing
  argument.
- **Carries.** All of Layer I (bound, block identity, ladder, slack), the
  relational-capacity coding theorem (genealogy has zero exponential tax in
  \(\mathbb H_\kappa^n\)), the Layer IIa realization theorems and remaining
  curvature-genericity conjecture, and \(n=2\) as embeddability. Presents the
  state equation as a *conditional* Layer IIb ideal, not a law.
- **Backed by.** The Lean development under `theory/lean/` (bound and block
  identity fully checked; Layer II algebra checked), the cited Bishop--Jones
  lower bound in `theory/RELATIONAL_CAPACITY_THEOREM.md`, and the
  packing/quartet mathematics of the spine.
- **Claim it makes about nature.** None that can fail — it is a converse plus an
  achievable form plus a realization geometry. Its risk is entirely
  mathematical (the curvature-genericity conjecture).

### Paper II — *Evolution as Active Geometry*

- **Domain.** Evolutionary biology, phylogenetics.
- **Carries.** The **Layer IIa** claim as its defensible spine — that
  phylogenetic data is hosted by hyperbolic tree geometry (seed-stable
  embeddings, curvature as a fixed design parameter, tree-defect near zero,
  \(n=2\) by embeddability) — and the **Layer IIb** state-equation/saturation
  claim as an *honestly ledgered open bet*, not a headline. Leads with the
  growth-class × tree-defect phase diagram, not with \(\kappa\).
- **Backed by.** The empirical pipelines under `validation/`, `experiments/`,
  `model/`, and `figures/`; the shared meter; the E-series protocol.
- **Claim it makes about nature.** Real and falsifiable, and currently mixed:
  host class comparatively supported but its decisive test (E9) unrun;
  saturation tested and failing its independent kill lines. Paper II must state
  this split rather than average over it.

A prospective **Paper III — *Language as Active Geometry*** (convergent
alphabets) would fill a non-biological substrate row; it is a pre-registered
future test, not part of the current artefact.

## The seams (joinery)

- **Kernel seam.** Paper II *cites* Paper I for the bound and block identity; it
  never re-derives them. The inequality is the license; the biology is an
  instance.
- **The genealogy seam.** Paper I proves that relational fidelity costs no
  exponential-order capacity in a real hyperbolic host under a local weighted
  process clock. Paper II therefore tests utilization and host class, not an
  assumed universal "price of genealogy." E2's equal-edge endpoint obstruction
  is a stronger synchronization issue, not host capacity.
- **The clock seam.** Paper I proves the radial concentration converse:
  near-capacity operation in a hyperbolic host forces almost all retained
  histories toward \(d(o,f(v))=c\tau(v)\), along a capacity-realizing sequence
  (or eventually when the growth limit exists), with exponent
  \(h\delta c-\Delta_{\rm cap}\). Here
  \(\Delta_{\rm cap}=ch-\beta=ch(1-\eta)\), not the dimensionless efficiency
  itself. Paper II's E6 tests this **conditional consequence of IIa+IIb**, not
  an unconditional claim that biological radius is information. A negative E6
  challenges the biological conjunction only after all theorem premises and
  unbiased clock-ball sampling are verified; independent host-class and
  utilization measurements identify which premise failed.
- **The \(n=2\) seam.** Paper I fixes \(n=2\) as a structural embeddability
  floor. Paper II therefore must **not** present \(n=2\) as an empirical
  discovery; a measured \(n>2\) (e.g. bacteria \(\approx 3.4\)) is a
  homogeneity-premise / reticulation scope statement (an E8 finding), not a
  refutation of the floor.
- **The \(\kappa\) seam.** Paper I defines \(\kappa\) as normalized (gauge-
  dependent) curvature and the state equation as conditional. Paper II may
  report a measured \(\kappa\) only with the isotropy premise flagged and E9's
  status attached; no absolute \(\kappa\) is a certified measurement today.
- **The instrument seam.** Both papers cite the same certified instruments
  (`tools/addressability_meter.py`: M2, M4, the Layer 0 growth-class gate in
  `tools/growth_class_gate.py`). The gate emits a class only when
  [`MEASURABILITY.md`](MEASURABILITY.md) marks the sample measurable. Neither
  paper may quote \(\eta\) to a precision the uncertified M3 magnitude
  estimator does not support (\(\pm20\%\) is not yet interpretable).
- **The finite-sample seam.** Paper I's rates are limsups. Paper II's
  matrices are finite. Layer 0 is the join: a matrix that fails the
  measurability predicate is not placed on the growth-class axis, and it
  cannot be read as evidence for or against saturation. Short codes are
  not approximate capacity-achieving codes.
- **The honesty seam.** [`CLAIMS.md`](CLAIMS.md) tags every claim and binds it
  to its artifact and layer; it is the contract that keeps Paper II's empirical
  risk out of Paper I.

## Self-verification

The artefact checks itself, and CI should run all five:

1. **Lean.** `cd theory/lean && lake build` — no `sorry`/`admit`;
   block-identity proofs use only `propext`, `Classical.choice`, `Quot.sound`.
2. **Meter.** `python3 -m unittest tests.test_addressability_meter` — 6 tests.
3. **Layer 0 gate.** `python3 -m unittest tests.test_growth_class_gate`.
4. **Doc/artifact registry.** `python3 tools/check_doc_artifacts.py` — fails if
   a tracked doc names a repository path not on disk (enforced on
   `theory`/`tools`/`tests`; `--all` advisory elsewhere).
5. **Claim ledger.** [`CLAIMS.md`](CLAIMS.md) is the single source of truth for
   what is THEOREM / IDENTITY / OPEN / CONVENTION / INSTRUMENT / EMPIRICAL.

## Status: near-finished, with the gaps named

- **Finished:** Layer I (proved, checked); the two-layer/three-sublayer
  framing; the instrument for the qualitative claims (M2, M4); the honesty
  apparatus.
- **Finished, Layer 0:** the growth-class measurability theorem, the Lean
  identities, and the instrument that refuses unmeasurable windows. The
  2×-span death of the regression gate is now a bound, not a surprise.
- **Finished, paper-level mathematics:** weighted/radial relational capacity of
  \(\mathbb H_\kappa^n\) equals its volume entropy (Skenderi lower bound plus
  the packing converse). The former unit-edge conjecture is refuted as stated.
- **Open, mathematics:** the curvature-genericity conjecture, plus the stronger
  equal-edge synchronization refinement for \(c\ge\varepsilon\) (not the
  definition of host capacity). Further Layer 0 finite-sample results
  (magnitude intervals, short-code signatures) are not opened here.
- **Open, instrument:** a certified M3 magnitude estimator and a runtime
  independence firewall.
- **Open, empirical — the highest-leverage gap:** E0 on real matrices, then
  E9. E9 remains the decisive IIa intervention; E0 is the filter that
  decides which matrices may be asked any growth-class question at all.

## Center of gravity — a stopping rule

A theorem of the form "if \(A\) and \(B\) then \(C\)" is a gift to the theory
and a debt to the experiment. Almost every recent Layer I/II theorem
lengthened that chain. The spine of Layers I and II remains **closed for
Paper I purposes**. Curvature genericity, the sector-gluing sharpness
lemma, the equal-edge refinement, and further asymptotic converses are
deferred.

Layer 0 is the exception, and the only one: a theorem whose hypotheses
*are* the measurement, and whose conclusion eliminates experiments rather
than adding premises. That exception is now discharged. Further Layer 0
work is admitted only if it likewise shortens the chain (a finite-sample
magnitude interval that retires M3-as-number; a short-code signature that
replaces asymptotic \(\eta\to 1\) as the IIb target). Another converse
that adds a clock, an offset, or an isotropy switch is not admitted.

The empirical queue, in order of immediacy:

1. **E0 — measurability audit, then the phase diagram.** Run the Layer 0
   predicate on real biological distance matrices. Place only the
   measurable subset on the growth-class × tree-defect plane. Unmeasurable
   matrices are recorded as such, not as failed classifications.
2. **E5 at small scale — runnable now.** Trained-hierarchy saturation with
   co-equal negative controls, on windows Layer 0 marks measurable.
3. **E4 rehearsal — public data.** Time-stamped serially sampled viruses as
   the pipeline rehearsal, labelled as such.
4. **E9 — the decisive test.** Needs GPU-scale embedding; highest IIa value.

The next unit of effort belongs to item 1.
