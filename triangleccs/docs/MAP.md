# The map on a frozen chart

Descent occupies the room; genomes encode it. TriangleCCS is how a query of
that encoding stays comparable. v10.9 (and successors) is how the encoding is
queried. This note is the long-horizon interface. It is not occupancy of the
top row, not E9, and not a filled atlas.

See `docs/ENCODER.md` for why the sextant, not a language model, is the
chart encoder. See `docs/VERSIONING.md` for epochs.

## Interface, not isomorphism

Registration is an \(O(2)\) conjugacy on the coarse polar backbone: tangent
map at the origin, a 2D backbone, two anchors (`triangleccs/datum/registration.py`).
The published object is a transform JSON (`schemas/transform.v1.json`):
backbone basis, tangent mean, `certified: false` until the freeze-gate
passes. Fitting a new atlas does not edit `Form`.

That mapping is not a metric isomorphism. Long-range JC saturates; the
encoder keeps block distinguishability. Leaf-contrastive training stretches
sibling genera. Extra latent dimensions are the map's resolution, not more
of the certified tree. The elegant correspondence is freeze-gate class
θ agreement on a stratified aligned panel, plus a named tip residual. If
after one global \(O(2)\) some phyla still need their own rotation, the
transform is a compromise, not a conjugacy, and the gate should refuse.

## Who runs what

A lab with a genome — long or short, seen or unseen — runs the encoder, then
the published transform, and receives an `Address`: θ, advisory r, Form
hash. The atlas **renders** species as those registered coordinates. A small
encoder is the compressed query; shrinking the model is in-scope for the
map.

The sextant is not that path. It consumes a named alignment, places JC
distances on the chart, fits or checks the transform, and measures
distortion. Users live in the atlas. The pin is a map reading on a frozen
chart, not a JC measurement of life.

Drop the map's radial head. Do not read live κ. Radius stays ADVISORY; κ
stays CONVENTION.

## Time horizon

Form epochs are rare. Changing κ, anchors, ε, or the radial-proxy name is a
new `Form.version` and a dual-address window, like a leap second.

Map retrains publish a new transform file. They never bump the Form.
Warm-start inheritance across atlas versions is `CIRCULAR` and cannot
satisfy the freeze-gate alone. The load-bearing witness is the sextant on a
shared aligned panel, not two nets that saw the same tree.

A free atlas and an improving atlas must share the Form hash. The improving
map may be better at tips. It may not have a private meridian.

## Scale

Millions of whole genomes are a map. The witness is a **stratified aligned
panel** — anchors, a backbone of representatives, dense local patches —
not dense JC on a large fraction of taxa. Most deep pairs on a biosphere
alignment are the balloon (v6 Test A): endpoints remain distinct, quartets
go quiet. Corroboration is coarse θ, not global distances. Extra dimensions
stay off the chart.

An operational probe of the live query path (v15.5 serving of the v10.9
map epoch) lives in `experiments/map_query/`. It records map behaviour
through Form tags; it is not a freeze-gate run and not a shared aligned
panel.

## Governance

The datum is a commons. The atlas may be corporate. Independent stewardship
after adoption is a horizon for the Form, not a result of Paper II. What a
steward would hold is the packing bound and the freeze-gate, not the
weights.

## What this does not claim

Saturation, a live curvature, LUCA at the origin, milliarcsecond precision,
or encoder θ as occupancy of exponential room. Those belong elsewhere, or
nowhere.
