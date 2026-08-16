# Versioning

Three layers change at different rates.

## Form epoch

`Form` is immutable for a given `version` string. Changing κ, dimension,
tokenizer, anchors, ε, or the radial-proxy *name* requires a new epoch and a
published transform from the previous epoch. The form hash is embedded in every
`Address`.

## Transform

A transform JSON maps a particular map realization (e.g. Atlas 129D) onto the
chart: backbone basis, tangent mean, anchor coordinates, `certified: false`
until the freeze-gate passes. Fitting a new atlas does **not** edit `triangleccs/datum/form.py`.
It publishes a new transform file.

Warm-start inheritance across atlas versions is tagged `CIRCULAR` and cannot
alone satisfy the freeze-gate. The load-bearing witness is an independent
sextant lineage (`triangleccs/sextant/place.py`: distances onto the chart, not
a genomic LM), when available.

## Atlas / map

The map (129D LM, Voronoi retrieval, UI) lives outside this repo. It *emits*
`Address` via registration. Extra dimensions beyond the 2D chart are operational
resolution of the map, not more of the certified tree. See `docs/ENCODER.md`.
