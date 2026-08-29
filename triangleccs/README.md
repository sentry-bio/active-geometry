# TriangleCCS

**A geodetic decoder datum for a tree-structured source on an overwriting tape.**

Not a globe of life. Not a language model. Not a second Lean spine.

Theorems live in [active-geometry](https://github.com/sentry-bio/active-geometry)
([throughline](https://github.com/sentry-bio/active-geometry/blob/main/theory/THROUGHLINE.md)).
Biosphere Atlas is a **consumer**: it registers onto this form. It does not
define the form.

## One sentence

Packing is the limit; genealogy is not packing; quartets ask whether the
received metric is still a tree; polar \(\mathbb H^2\) is the candidate chart
in which those two uses of capacity can coincide. This package freezes that
chart as a WGS84-like datum. The encoder that belongs here is a **sextant**
(distances onto the chart), not a genomic language model — see
[`docs/ENCODER.md`](docs/ENCODER.md). How a map (v10.9) registers over a
long horizon is [`docs/MAP.md`](docs/MAP.md).

## Install

```bash
pip install -e ".[dev]"
pytest
python examples/register_atlas.py
python examples/run_balloon.py
python experiments/map_query/run_benchmark.py
python tools/check_doc_artifacts.py
```

## Public API

```python
from triangleccs import Form, make_address, place_sequences, DEFAULT_FORM

form = Form()  # κ=1.25 CONVENTION, dim=2 inhabit H², r advisory
print(form.summary())
```

## What is not claimed

- κ = (h ln 2)² is not a theorem of this datum; κ is a frozen **CONVENTION**.
- Radius is not accumulated information; it is an **ADVISORY** named proxy.
- Chart origin is not LUCA.
- Cross-instrument θ agreement is reproducibility within the imposed model,
  not information-tightness or saturation.
- Curvature genericity and the state equation are overlays (empty in v1).

## Layout

See `CONSTITUTION.md` and `CLAIMS.md`. Four pieces map to
`triangleccs/packing/`, `triangleccs/tape/`, `triangleccs/classifier/`,
`triangleccs/chart/` + `triangleccs/datum/`. The sextant lives in
`triangleccs/sextant/`. The shared metric is `triangleccs/metric.py`.

## License

MIT.
