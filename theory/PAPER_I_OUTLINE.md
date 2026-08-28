# Paper I outline — *The Addressability Limit*

Draft for the private manuscript that replaces the Zenodo record
*A Geometric State Equation for Information-Generating Hierarchies*
(concept DOI [10.5281/zenodo.19381557](https://doi.org/10.5281/zenodo.19381557),
v3 [10.5281/zenodo.21383132](https://doi.org/10.5281/zenodo.21383132))
before arXiv. Not a TeX file. The proof body is
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md); this page is the
publication contract.

Paper II cites this paper and instantiates it. It does not re-derive it
([`PAPER_II_OUTLINE.md`](PAPER_II_OUTLINE.md)).

## Title

*The Addressability Limit: A Packing Bound for Information-Generating
Hierarchies*

The load-bearing object is a packing inequality in a proper metric space.
The state equation is a detachable cost, named in a remark, not in the
title. That is the replacement: the earlier public face was the equality;
this manuscript is the ceiling that equality was pretending to be.

## Contribution

A process that retains distinguishable histories at finite radial cost
cannot outrun the packing entropy of its host. The same count, with the
quantifiers flipped, is a trichotomy: exponential room, exponentially
growing addresses, or forgetting. Block addressability is exactly packing.
Preserving genealogy is a strictly stronger task; in a weighted hyperbolic
host the exponential-order tax on that task is zero. Quartets classify
tree-ness, orthogonally. Polar \(\mathbb H^2\) is the minimal-dimensional
candidate in a stated embedding class, given an isotropy axiom that is
named as a gauge.

No biology in the load-bearing argument. No claim that nature saturates.
No claim that a bake-off forces the host.

## arXiv category

**Primary: `math.MG` (Metric Geometry).** The theorem is a packing bound
for \(\varepsilon\)-separated sets in proper metric spaces, with packing
entropy as the host invariant. Volume comparison (Bishop–Gromov at fixed
dimension) is a corollary engine, not a Riemannian premise of the ceiling.

**Secondary: `cs.IT` (Information Theory).** The operational reading is
block address capacity, then a ladder of constrained capacities. The
audience overlap is Shannon-style converses, not a channel-coding paper
that assumes a named discrete channel.

**Do not submit as**

- `q-bio.*` — that is Paper II;
- `cs.LG` — hyperbolic embeddings as a model class are the literature this
  paper is *not* joining;
- `cond-mat.stat-mech` — the retired title's "state equation" neighborhood;
- `math.DG` as primary — the limit is curvature-free. Cross-list `math.DG`
  only if the Heintze/A3 realization section stays long in the submitted
  TeX.

MSC (if asked): 51K05 (distance geometry), 53C23 (global geometric
methods; hyperbolicity), 94A24 (coding theorems).

This repo is the private source until the TeX is cut. On posting: arXiv
first under those categories, then a new version of the Zenodo record
under this title, so the concept DOI does not keep advertising a forced
equality.

## What this manuscript keeps from the earlier record

- The opening tension: remembering while creating has nowhere polynomial
  to go at finite address cost.
- Finite rate as a *condition* of the ceiling, not a standing property of
  every representation.
- Compact exhaustion: every finite-depth census lives on a compact ball;
  the rate statement is the \(R\to\infty\) limit, remainder \(O(1/R)\).
- Derivation versus demonstration as a firewall. The result is a candidate
  principle about systems that satisfy the premises.
- Isotropy as a gauge. Anisotropy is a band, not an escape.
- Lean scoped to what it actually checks.

## What it refuses to carry forward

- Boxing \(h_{\mathrm{vol}}=h_{\mathrm{eff}}\ln 2\) as forced by the same
  axioms as the packing inequality. That equality is losslessness binding
  in a cost, Layer IIb, [`IIB_CONTRACT.md`](IIB_CONTRACT.md).
- A2 (unimodal channel, interior \(h^\ast\)) as load-bearing for the
  limit. The channel clock is not packing.
- Fisher–Rao location-scale geometry as a second derivation. Depth versus
  one-step (sphere versus \(\mathbb H\)) may appear as a remark; it does
  not prove the bound.
- Title-level \(\kappa\), or \(n\) back-solved from curvature.
- An experiment that "forces hyperbolic geometry." Corollary 4.3 already
  excludes polynomial hosts at finite rate and fixed dimension. E9 is a
  finite-sample shadow of that corollary at matched packing, not a pillar
  of this paper ([`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md)).

## How the paper reads

One argument. The four throughline pieces
([`THROUGHLINE.md`](THROUGHLINE.md)), in order.

- **§1** — remembering while creating; the trichotomy; the sharp outside
  (processes that forget are not governed, and that is the content).
- **§2** — the packing bound \(\beta\le c\,h_{\mathrm{pack}}\); the block
  identity; compact exhaustion; \(\varepsilon\)-robustness of the
  exponential *rate*.
- **§3** — the ladder; the balloon; weighted \(\mathbb H_\kappa^n\) has
  zero exponential-order relational tax (Theorem 4.4).
- **§4** — quartets; undecidable \(\neq\) confidently wrong; \(\delta\)
  orthogonal to packing.
- **§5** — chart: \(n=2\) as embeddability floor; A3 \(\Rightarrow\) real
  \(\mathbb H^{d+1}\) inside Heintze; polar radius/angle as encoder
  choices. Isotropy named as gauge. Uniform mixing raises dimension;
  structured mixing breaks isotropy.
- **§6** — firewall. What is derived; what would count as a demonstration;
  the state equation as a conditional ideal, not a law. No empirical
  tables.

Lean inventory and the Bishop–Jones citation for Theorem 4.4 belong in
appendices, with the same honesty as the earlier record: algebra is
checked; the geometry lemmas are classical and cited.

## Joints with Paper II (named, not used)

Paper II needs four licenses, already on the throughline: room, names
versus paths, quiet versus wrong, illegal identification. This paper
supplies them. It does not place a clade, fit a \(\kappa\), or run E9.

## Status

Mathematics: the packing bound and finite-block identity are proved (the
convergent-rate bound and finite-radius identity are Lean-checked). The
capacity ladder is a hierarchy of definitions; its rung values remain open
except for block capacity and the host-specific weighted relational
supremum of Theorem 4.4, whose lower bound is paper-level and cited.
Theorem 7.1 remains a paper sketch conditional on A3. Cutting TeX requires
no stronger claim than this boundary.
