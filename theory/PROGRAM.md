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
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md); the tests in
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
  Layer IIa realization theorems and the two open conjectures (relational
  capacity, curvature genericity), and \(n=2\) as embeddability. Presents the
  state equation as a *conditional* Layer IIb ideal, not a law.
- **Backed by.** The Lean development under `theory/lean/` (bound and block
  identity fully checked; Layer II algebra checked; conjectures stated, not
  proved) and the packing/quartet mathematics of the spine.
- **Claim it makes about nature.** None that can fail — it is a converse plus an
  achievable form plus a realization geometry. Its risk is entirely
  mathematical (the two conjectures).

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
- **Open, mathematics:** the relational-capacity conjecture (subcritical
  achievability; E2 found only an endpoint obstruction) and the
  curvature-genericity conjecture.
- **Open, instrument:** a certified M3 magnitude estimator and a runtime
  independence firewall.
- **Open, empirical — the highest-leverage gap:** E9, the matched-capacity
  Euclidean-vs-hyperbolic realization test. The biological claim is Layer IIa,
  and E9 is the only direct IIa intervention; it is designed and unrun.
