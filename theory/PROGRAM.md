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
  \(C_{\mathrm{block}}=c\,h_{\mathrm{pack}}\); the constrained-capacity
  ladder as nested admissibility classes; the slack decomposition; the
  trichotomy. Proven as packing geometry. Mentions no curvature.

- **Host-class closure of the relational rung.** Theorem 4.4: in
  \(\mathbb H_\kappa^n\), under a local weighted process clock, the
  supremum relational capacity equals block capacity. This is a hyperbolic
  host theorem, not a Layer I universal result. It is stated next to the
  ladder because that is where the question arises.

- **Layer IIa — host class (realization).** *Which* geometry hosts the data:
  hyperbolic, at what dimension, given named symmetry. The space-form argument
  (given isotropy), Theorem 7.1 (Heintze isotropy given axiom A3), \(n=2\) as
  an embeddability floor. The **better-supported** biological claim is
  *occupancy* of exponential tree-like room, not a bake-off that forces the
  host. A3 is asserted, not measured. E9 is Corollary 4.3's finite-sample
  shadow at matched packing.

- **Layer IIb — saturation.** *Whether* a process fills its budget, giving the
  state-equation equality. The **harder, less-supported** claim; every
  independent test of the *retired plot* (sequence-\(\kappa\) vs tree-derived
  \(h\)) sits below its kill line. The well-posed bet is
  [`IIB_CONTRACT.md`](IIB_CONTRACT.md): \(\eta\) on the process after
  exponential occupancy, mixed cell out. It is a sequel, not a chapter of
  Paper II.

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

This is the private manuscript that replaces the Zenodo record *A Geometric
State Equation for Information-Generating Hierarchies* before arXiv. The
publication contract is [`PAPER_I_OUTLINE.md`](PAPER_I_OUTLINE.md).

- **Domain.** Metric geometry, information theory. No biology in the load-bearing
  argument.
- **arXiv.** Primary `math.MG`; secondary `cs.IT`. Not `q-bio.*`, not
  `cs.LG`, not `cond-mat.stat-mech`. Cross-list `math.DG` only if the
  Heintze section stays long in the submitted TeX.
- **Carries.** All of Layer I (bound, trichotomy, block identity, ladder
  as nested classes, slack), the relational-capacity coding theorem as a
  *host-class* identity (genealogy has zero exponential-order tax in
  \(\mathbb H_\kappa^n\) under a weighted clock; supremum, not
  attainment), the Layer IIa realization theorems including Theorem 7.1
  (Heintze isotropy, conditional on axiom A3), and \(n=2\) as
  embeddability. Presents the state equation as a *conditional* Layer IIb
  face, not a law, and not in the title.
- **Backed by.** The Lean development under `theory/lean/` (bound and block
  identity fully checked; Layer II algebra checked), the cited Bishop--Jones
  lower bound in `theory/RELATIONAL_CAPACITY_THEOREM.md`, and the
  packing/quartet mathematics of the spine.
- **Claim it makes about nature.** None that can fail — it is a converse plus an
  achievable form plus a realization geometry. Its remaining mathematical
  risk is the paper sketch of Theorem 7.1 (not Lean-checked) and the
  equal-edge synchronization refinement. Axiom A3 is not a Paper I claim
  about nature. E9 is not a Paper I pillar.

### Paper II — *Evolution as Active Geometry*

- **Domain.** Evolutionary biology, phylogenetics. The living instance of
  Layer IIa (host class) and the open Layer IIb bet (saturation). It does
  not carry a coordinate-system atlas, a self-calibrating information clock,
  or a claim that selection fills packing capacity.
- **Carries.** The **Layer IIa** claim as occupancy of the program's figure —
  growth class \(\times\) tree defect, on a representation that is not the
  inferred tree — and the **Layer IIb** state-equation/saturation claim as
  an *honestly ledgered open bet*, not a headline. The working outline is
  [`PAPER_II_OUTLINE.md`](PAPER_II_OUTLINE.md): descent occupies exponential
  room; genomes encode it; sequence metrics lose the tree by erasure or by
  reticulation; sequence-\(\kappa\) was a plot of the encoding.   A3 and saturation unrun. Occupancy, not force: Corollary 4.3 already
  excludes the bottom row at finite rate; E9 is an optional shadow.
- **Backed by.** The empirical pipelines under `validation/`, `experiments/`,
  `model/`, and `figures/`; the shared meter; the E-series protocol.
- **Claim it makes about nature.** Real and falsifiable, and currently mixed:
  the process occupies the top row of the figure on independent functional
  geometry; mixing scars are
  two independent HGT instruments, with Test B still the gate on calling a
  real-data floor mixing; saturation tested and failing its independent kill
  lines. In the intended hyperbolic host those failures are utilization, not
  relational tax (Theorem 4.4). Paper II must state this split rather than
  average over it. A matched-parameter Euclidean win at extra \(d\) does not
  kill occupancy.

A prospective **Paper III — *Language as Active Geometry*** (convergent
alphabets) would fill a non-biological substrate row; it is a pre-registered
future test, not part of the current artefact.

## The seams (joinery)

- **Kernel seam.** Paper II *cites* Paper I for the bound, the trichotomy,
  and the block identity; it never re-derives them. The inequality is the
  license; the biology is an instance. The manuscript contract is
  [`PAPER_I_OUTLINE.md`](PAPER_I_OUTLINE.md).
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
  report a measured \(\kappa\) only with the isotropy premise flagged; no
  absolute \(\kappa\) is a certified measurement today.
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
  supplied by construction — not a host-class selection test (occupancy
  of the figure is; E9 is a packing-matched illustration). The two
  training axes — NCBI ranks for quartets, genome size for
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
- **Open, mathematics:** the equal-edge synchronization refinement for
  \(c\ge\varepsilon\) (not the definition of host capacity). Theorem 7.1 is a
  paper sketch, not Lean. Axiom A3 is empirical (L2a.6), not a remaining
  geometric conjecture.
- **Open, instrument:** a certified M3 magnitude estimator and a runtime
  independence firewall. The growth-class gate refuses short radial
  windows; that is a property of the instrument, recorded in
  [`MEASURABILITY.md`](MEASURABILITY.md).
- **Open, empirical — occupancy, not force:** Paper II's figure on an
  independent representation (growth class \(\times\) tree defect). The
  biological claim is Layer IIa occupancy. Corollary 4.3 already excludes
  the polynomial row at finite rate; E9 is optional.

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
   \(\kappa\). This is Paper II's figure
   ([`PAPER_II_OUTLINE.md`](PAPER_II_OUTLINE.md)). A first pass on
   tree-patristic matrices refused 11/24 windows and made \(\delta\)
   vacuous by construction; that is an instrument result, not a placement
   of biology.
2. **E5 at small scale — stopped on a design flaw.** Trained-hierarchy
   saturation remains the near-term IIb item of
   [`IIB_CONTRACT.md`](IIB_CONTRACT.md), but the current embedder is not
   an exponential-capacity chart of the generator. Report:
   [`experiments/E5_DESIGN_FLAW.md`](../experiments/E5_DESIGN_FLAW.md).
   Do not report \(\eta\). Redesign before rerunning.
3. **E4 rehearsal — public data.** Time-stamped serially sampled viruses as
   the pipeline rehearsal, labelled as such.
4. **E9 — optional illustration.** Matched-*capacity* Euclidean vs
   hyperbolic, never matched parameter count. Finite-sample shadow of
   Corollary 4.3; not the test that forces the host.

The next unit of effort belongs to item 1, not to another converse, and
not to a new domain. Cutting Paper I TeX is a manuscript act
([`PAPER_I_OUTLINE.md`](PAPER_I_OUTLINE.md)), not a new theorem.
