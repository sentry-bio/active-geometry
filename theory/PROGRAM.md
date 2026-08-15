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

## Three layers

The program is stratified so that a failure in one layer cannot contaminate
another. The definitions live in [`CLAIMS.md`](CLAIMS.md); the mathematics in
[`ADDRESSABILITY_KERNEL.md`](ADDRESSABILITY_KERNEL.md) and
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md); the relational coding theorem
in [`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md); the tests in
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md).

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
  (`tools/addressability_meter.py`: M2, M4, the growth-class gate). Neither may
  quote \(\eta\) to a precision the uncertified M3 magnitude estimator does not
  support (\(\pm20\%\) is not yet interpretable).
- **The honesty seam.** [`CLAIMS.md`](CLAIMS.md) tags every claim and binds it
  to its artifact and layer; it is the contract that keeps Paper II's empirical
  risk out of Paper I.

## Self-verification

The artefact checks itself, and CI should run all four:

1. **Lean.** `cd theory/lean && lake build` — 3017 jobs, no `sorry`/`admit`;
   block-identity proofs use only `propext`, `Classical.choice`, `Quot.sound`.
2. **Meter.** `python3 -m unittest tests.test_addressability_meter` — 6 tests.
3. **Doc/artifact registry.** `python3 tools/check_doc_artifacts.py` — fails if
   a tracked doc names a repository path not on disk (enforced on
   `theory`/`tools`/`tests`; `--all` advisory elsewhere).
4. **Claim ledger.** [`CLAIMS.md`](CLAIMS.md) is the single source of truth for
   what is THEOREM / IDENTITY / OPEN / CONVENTION / INSTRUMENT / EMPIRICAL.

## Status: near-finished, with the gaps named

- **Finished:** Layer I (proved, checked); the two-layer/three-sublayer
  framing; the instrument for the qualitative claims (M2, M4, growth-class
  gate); the honesty apparatus.
- **Finished, paper-level mathematics:** weighted/radial relational capacity of
  \(\mathbb H_\kappa^n\) equals its volume entropy (Skenderi lower bound plus
  the packing converse). The former unit-edge conjecture is refuted as stated.
- **Open, mathematics:** the curvature-genericity conjecture, plus the stronger
  equal-edge synchronization refinement for \(c\ge\varepsilon\) (not the
  definition of host capacity).
- **Open, instrument:** a certified M3 magnitude estimator and a runtime
  independence firewall.
- **Open, empirical — the highest-leverage gap:** E9, the matched-capacity
  Euclidean-vs-hyperbolic realization test. The biological claim is Layer IIa,
  and E9 is the only direct IIa intervention; it is designed and unrun.

## Center of gravity — a stopping rule

Effort in this program flows downhill toward mathematics, because theorems can
be produced on demand and experiments cannot. Every recent theorem also
lengthened the dependency chain from proof to biological meaning (E6 now
requires verified host class, measured utilization, and a calibrated clock
before it means anything). That asymmetry is the architecture working — Paper I
may advance while Paper II waits — but unchecked it is the failure mode of
mathematically beautiful programs that lose contact with their subject.

Accordingly, the spine is declared **closed for Paper I purposes**. The
remaining mathematical items (curvature genericity, the sector-gluing
sharpness lemma, the equal-edge refinement, Lean formalization of asymptotic
corollaries) are not load-bearing for any biological claim and are deferred.

**No new theorems until a measurement has run.** The empirical queue, in order
of immediacy:

1. **The phase diagram on real matrices — runnable now.** The growth-class
   gate and the δ-meter are the two *certified* instruments; they need no
   \(\eta\) precision and no \(\kappa\). Place real biological distance
   matrices on the growth-class × tree-defect plane.
2. **E5 at small scale — runnable now.** Trained-hierarchy saturation with
   co-equal negative controls.
3. **E4 rehearsal — public data.** Time-stamped serially sampled viruses as
   the pipeline rehearsal, labelled as such.
4. **E9 — the decisive test.** Needs GPU-scale embedding; highest value.

The next unit of effort belongs to item 1, not to another converse.
