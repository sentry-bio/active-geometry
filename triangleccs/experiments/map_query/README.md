# Live map query probe (v10.9 epoch / v15.5 serving)

Operational reading of the Biosphere Atlas query path **through TriangleCCS**.
It does not certify Form and it does not treat the globe as occupancy.

## What is harnessed today

```
DNA ──► POST /identify ──► V15.5 encode (129D Poincaré)
                         ──► geodesic hierarchical placement (Domain→Species)
                         ──► conformal zone (accept / escalate / fallback)
                         ──► atlas_r, atlas_theta, 3D PCA coords, live κ

DNA ──► POST /predict   ──► domain/family/genus heads
                         ──► 3D PCA (free) + 129D tangent (API key)

UI (biosphereatlas.com) ──► IDENTIFY calls /identify, then looks the
                            returned species up in ball_data.json (v9).
```

The live API reports **v15.5** (OpenAPI 15.5.5). TriangleCCS docs name the
map epoch **v10.9**. Same role: 129D LM + Voronoi/geodesic retrieval + UI.
The static ball is still a **v9 consensus** index (121,351 genomes, display
κ = 1.2453). Live κ on the GPU was 1.2369 when this probe was designed.

## What is not harnessed

- No `Form` hash on responses.
- No published 129D → polar-chart transform (`examples/atlas_transform.v1.json`
  is an 8D fixture).
- No sextant: there is no shared *aligned* panel.
- Public path still returns `atlas_r` and live κ. MAP.md says drop the
  radial head and do not read live κ.
- `/identify` does not emit a TriangleCCS `Address`.

## What this probe asks

| Probe | Question | TriangleCCS reading |
|---|---|---|
| Domain panel | Do 16S/18S from three domains land in the right domain? | Map quality, not occupancy |
| `atlas_theta` | Does the map's private polar angle cluster by domain? | Private meridian until registration |
| 129D tangent | After dropping a radial head, does SVD 2D θ separate domains? | CANDIDATE θ, Form κ CONVENTION |
| Length ladder | Same E. coli 16S at 100…full bp | θ should be stabler than r; r is ADVISORY |
| Siblings | Enterobacteriaceae vs distant bacteria | Leaf-contrastive over-separation / confusion |
| Nulls | shuffled / random / poly-AT / short | OOD should escalate, not look like a species |
| Mash k-mers | Second map vs Poincaré of the 129D vectors | Independent string metric; not a sextant |
| Freeze-gate | Can we flip θ to certified? | No — no aligned sextant witness |

## Run

```bash
cd triangleccs
python experiments/map_query/run_benchmark.py
```

Optional: `BIOSPHERE_API_KEY` (otherwise the public query-panel key is
discovered from `query-panel.js`). The accession panel is fetched from NCBI
on first run and cached as ignored `sequences.fasta`. Results are generated
locally at `results/latest.json` and are not committed.

## Honesty

κ on Form stays 1.25 CONVENTION. Chart origin is not LUCA. Taxonomy
accuracy is not saturation. 3D PCA is a decoder view of the map, not the
polar chart. Mash distances are an INSTRUMENT, not a topology supervisor.
