# The encoder that belongs to this datum

A genomic language model is **not** the optimal tool for an MDL expression of
the diversity of life. It is a useful **map** of strings. This repository's
encoder is a **sextant**: pairwise distances onto a frozen polar chart.

## Two description lengths

**MDL of sequences.** Compress genomes. A language model, a tree-aware
compressor, or gzip can all be candidates. Extra latent dimensions help
reconstruct tokens (the MLM packing problem). That objective is not a tree.

**MDL of diversity.** Compress *who is related to whom, and who is a distinct
address*. The short description is a tree (or its quartets) plus an
ε-packing of places. Phylogenetic likelihood *is* that MDL under a
substitution model. Packing is the MDL of distinguishable endpoints.

Those are different codes. A contrastive genomic encoder is a lossy hash of
strings that often *correlates* with the tree. Correlation is not the code.

## Why the globe reading failed

Active Geometry's throughline already splits the pieces a language model
tried to glue together:

- InfoNCE temperature and κ are non-identifiable, so a trained curvature is
  not a measurement.
- HEX (tree) and MLM (reconstruction) plateau in different latent dimensions.
  The embeddability floor of a tree is \(n=2\); extra dimensions are encoder
  slack and sequence reconstruction, not more of the certified tree.
- The balloon: sequence distance can keep endpoints distinct after quartets
  stop being readable. The metric a genomic encoder sees is the tape, not
  the genealogy.
- Taxonomy-supervised InfoNCE still trains against a tree of life. It is not
  an independent instrument.

So the Atlas globe is a **decoder view of a map**, not the MDL of diversity
and not a filled \(\mathbb H^2\).

## What is optimal, given that

| Role | Tool | Status |
|---|---|---|
| Observations | aligned or unaligned sequences | data |
| Channel | substitution / overwriting tape (balloon) | INSTRUMENT |
| Relations | quartets / additive tree | THEOREM (four-point) |
| Place | packing in a pointed metric | THEOREM (cite AG) |
| Chart | polar \(\mathbb H^2\), κ CONVENTION, r ADVISORY | this Form |
| Map / registration | genomic LM, Mash, k-mers, … | consumer of Form |
| Sextant | distances → chart, independent of the LM | this package |

The genomic encoder remains the right *operational* tool for placing
unaligned or metagenomic strings onto the chart **once the form exists**.
It is the wrong tool to *define* the form, to set κ, to certify radius, or
to stand in for the tree.

The MDL-relevant encoder is therefore the sextant in
`triangleccs/sextant/place.py`: JC (or any named) distances into polar
\(\mathbb H^2\), with δ, resolvability, block-separation, and residual filled
on every `Address`. A 129-dimensional language model may register onto that
chart. It must not author it.

## What the sextant reads

The sextant reads sequences. It does not read them the way a genomic encoder
does.

It consumes **aligned** bases, turns them into a named pairwise metric
(Hamming / JC69 in v0), and places that metric on the frozen polar chart by
hyperbolic law of cosines. That is a reading: the sequences move the needle.
It is a **compass reading**, not a survey of the coastline.

It will not, and should not:

- embed unaligned windows, short reads, or metagenomic queries;
- serve 10⁵ genomes as a nearest-neighbour index;
- invent angles that a contrastive loss prefers to the distances;
- author κ, n, or radius.

Those jobs belong to the **map** (v10.9 and kin). The map registers onto the
chart; the sextant is the non-neural reference that measures the map's
distortion against distance-faithful placement. The long-horizon interface
(epochs, query path, scale, governance) is `docs/MAP.md`. Two nets trained on the same
tree agreeing was always a soft form of independence. A JC sextant cannot
share a network's biases because it is not a network.

The gating requirement for that measurement is a **shared reference set**:
marker-gene alignments the sextant can place, tokenized windows the map can
embed, overlapping on real genomes. Until that set exists, a v10.9 registration transform cannot be fit or
conformance-checked. That gap is operational, not a reason to teach the
sextant to "read" like the map.

Fine-scale disagreement is expected and is the point. A leaf-contrastive map
over-separates sibling genera; the sextant preserves JC distances. Coarse
backbone θ should agree (conformance gate). Tip residual is the map's
characterized distortion, now measured against a distance-faithful ruler
instead of another net. Drop the map's radial head and ignore its live κ:
radius stays advisory on the chart; κ stays CONVENTION on the Form.

CONSTITUTION clause 2 still applies: the sequence metric places points; it
does not supervise topology. Quartets ask whether the received metric is
still a tree. The balloon says that at depth it often is not. So the sextant
is not "the true tree reader" either. It is the incorruptible *placement* of
a named tape metric onto the datum.

## What this does not claim

- That neighbor-joining or JC is the true process.
- That the sextant saturates host capacity.
- That laboratories should stop training genomic models. They should stop
  treating those models as the globe, the curvature, or the MDL of life.
