# Active Geometry — the program, whole

This is the administrative map of the whole artefact: the layered theory, the
two manuscripts it supports, the seams that join them, and the machinery that
makes the whole self-verifying. The shorter conceptual front door is
[`THROUGHLINE.md`](THROUGHLINE.md): limit, ladder, quartet classifier, chart.
Each part below points to the document that carries the detail.

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
  hyperbolic, at what dimension, and is that *forced*. The space-form argument
  (given isotropy), the curvature-genericity conjecture now narrowed to
  isotropy as exchangeability of the generator, \(n=2\) as an embeddability
  floor. The **better-supported** biological claim — but its premise,
  isotropy, is asserted, not measured. Decisive test: E9.

- **Layer IIb — saturation.** *Whether* a process fills its budget, giving the
  state-equation equality. The **harder, less-supported** claim; every
  independent test to date sits below its kill line.

The four-point (Buneman ≡ Gromov) classifier is an independent bridge to IIa:
it decides tree-ness and measures reticulation, but calibrates nothing.

These are the only domains. Finite-sample refusal of the growth-class gate
is an instrument limit ([`MEASURABILITY.md`](MEASURABILITY.md)). Radial
concentration is a corollary of IIa+IIb. Neither is a layer, a manuscript,
or a reason to enlarge the map.

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

- **Domain.** Evolutionary biology, phylogenetics. The living instance of
  Layer IIa (host class) and the open Layer IIb bet (saturation). It does
  not carry a coordinate-system atlas, a self-calibrating information clock,
  or a claim that selection fills packing capacity.
- **Carries.** The **Layer IIa** claim as its defensible spine — that a
  *representation* of phylogenetic data (not the inferred tree read as its
  own metric) is hosted by hyperbolic tree geometry — and the **Layer IIb**
  state-equation/saturation claim as an *honestly ledgered open bet*, not a
  headline. The honest figure is the growth-class × tree-defect phase
  diagram on an independent representation, not \(\kappa\).
- **Backed by.** The empirical pipelines under `validation/`, `experiments/`,
  `model/`, and `figures/`; the shared meter; the E-series protocol.
- **Claim it makes about nature.** Real and falsifiable, and currently mixed:
  host class comparatively supported but its decisive test (E9) unrun;
  saturation tested and failing its independent kill lines. In the intended
  hyperbolic host those failures are utilization, not relational tax
  (Theorem 4.4). Paper II must state this split rather than average over it.

A prospective **Paper III — *Language as Active Geometry*** (convergent
alphabets) would fill a non-biological substrate row; it is a pre-registered
future test, not part of the current artefact.

## The seams (joinery)

- **Kernel seam.** Paper II *cites* Paper I for the bound and block identity; it
  never re-derives them. The inequality is the license; the biology is an
  instance.
- **The genealogy seam.** Paper I proves that the *supremum* relational
  capacity of a real hyperbolic host under a local weighted process clock
  equals block capacity. Paper II therefore tests utilization and host
  class, not an assumed universal "price of genealogy." A measured
  \(\eta_{\mathrm{block}}<1\) is utilization slack, not tax, **only under
  a hyperbolic-host assumption**. On a finite-alphabet sequence metric at
  depth the same gap can be relational tax. E2's equal-edge endpoint
  obstruction is a stronger synchronization issue, not host capacity.
- **The clock seam.** Radial concentration is a *corollary* of IIa+IIb, not
  a load-bearing claim. If a hyperbolic host is near capacity, histories
  concentrate on the outer clock shell. E6 tests that consequence only after
  host class and utilization are independently known. It does not define
  the theory and it does not license a biological radius-as-information
  claim on its own.
- **The \(n=2\) seam.** Paper I fixes \(n=2\) as a structural embeddability
  floor. Paper II therefore must **not** present \(n=2\) as an empirical
  discovery; a measured \(n>2\) (e.g. bacteria \(\approx 3.4\)) is a
  homogeneity-premise / reticulation scope statement (an E8 finding), not a
  refutation of the floor.
- **The \(\kappa\) seam.** Paper I defines \(\kappa\) as normalized (gauge-
  dependent) curvature and the state equation as conditional. Paper II may
  report a measured \(\kappa\) only with the isotropy premise flagged and E9's
  status attached; no absolute \(\kappa\) is a certified measurement today.
- **The CCS instrument seam.** The Poincaré encoder remnant
  (`experiments/minimal_encoder/`, `model/`) is a Layer IIa *instrument*,
  not a Layer IIb measurement and not Paper II's claim. Polar
  \(\mathbb H^2\) is the embeddability floor; interpreting radius as
  depth and angle as divergence is an encoder modeling choice, not a
  capacity-theorem coordinate system. Frozen \(\kappa\) is the right
  response to InfoNCE degeneracy, while \(\kappa=5/4\) from genetic-code
  entropy is not a theorem; taxonomy quartets train topology and do not
  calibrate curvature. Seed-stable Procrustes (unique up to \(O(2)\),
  unless orientation is fixed) is reproducibility **within the imposed
  model** — \(\mathbb H^2\), curvature, and both training axes are
  supplied by construction — not a host-class selection test (that is
  E9). The two training axes — NCBI ranks for quartets, genome size for
  radius — are a primitive independence split, better than
  patristic-for-both, still not a representation metric independent of
  the tree of life. Radius-as-information, a filled atlas, and a
  self-calibrating clock remain IIb and stay out.
- **The instrument seam.** Both papers cite the same certified instruments
  (`tools/addressability_meter.py`: M2, M4, the growth-class gate). The gate
  refuses short radial windows; that refusal is an instrument limit
  ([`MEASURABILITY.md`](MEASURABILITY.md)), not a fourth layer of the
  theory. Neither paper may quote \(\eta\) to a precision the uncertified
  M3 magnitude estimator does not support (\(\pm20\%\) is not yet
  interpretable).
- **The honesty seam.** [`CLAIMS.md`](CLAIMS.md) tags every claim and binds it
  to its artifact and layer; it is the contract that keeps Paper II's empirical
  risk out of Paper I.

## Self-verification

The artefact checks itself, and CI should run all four:

1. **Lean.** `cd theory/lean && lake build` — no `sorry`/`admit`;
   block-identity proofs use only `propext`, `Classical.choice`, `Quot.sound`.
2. **Meter.** `python3 -m unittest tests.test_addressability_meter tests.test_growth_class_gate`.
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
- **Open, mathematics:** the curvature-genericity conjecture, narrowed to
  whether isotropy follows from exchangeable branching (space-form + packing
  already force \(\mathbb H_\kappa^n\) given isotropy), plus the stronger
  equal-edge synchronization refinement for \(c\ge\varepsilon\) (not the
  definition of host capacity).
- **Open, instrument:** a certified M3 magnitude estimator and a runtime
  independence firewall. The growth-class gate refuses short radial
  windows; that is a property of the instrument, recorded in
  [`MEASURABILITY.md`](MEASURABILITY.md).
- **Open, empirical — the highest-leverage gap:** E9, the matched-capacity
  Euclidean-vs-hyperbolic realization test. The biological claim is Layer IIa,
  and E9 is the only direct IIa intervention; it is designed and unrun.

## Center of gravity — a stopping rule

The theory is the four-piece dependency graph in
[`THROUGHLINE.md`](THROUGHLINE.md): limit, constrained-capacity ladder,
quartet classifier, and polar hyperbolic chart. The balloon is the separation
between room and room that genealogy fits; the coordinate-system construction
is the candidate host where Theorem 4.4 closes that separation. That spine is
**closed**. Radial concentration is a corollary of IIa+IIb. The growth-class
gate's short-window refusal is an instrument limit. Neither is a new
foundation, and neither is a reason to keep proving.

**No new theorems until a measurement has run.** The empirical queue stays
inside Paper II's domain. In order of immediacy:

1. **The phase diagram on an independent representation.** Growth class
   (refusing windows the gate cannot speak on) and tree defect, from a
   metric that is not the inferred tree. No \(\eta\) precision, no
   \(\kappa\). A first pass on tree-patristic matrices refused 11/24
   windows and made \(\delta\) vacuous by construction; that is an
   instrument result, not a placement of biology.
2. **E5 at small scale — runnable now.** Trained-hierarchy saturation with
   co-equal negative controls (Layer IIb, labelled as such).
3. **E4 rehearsal — public data.** Time-stamped serially sampled viruses as
   the pipeline rehearsal, labelled as such.
4. **E9 — the decisive test.** Matched-capacity Euclidean vs hyperbolic.
   Needs GPU-scale embedding; highest IIa value.

The next unit of effort belongs to item 1, not to another converse, and
not to a new domain.
