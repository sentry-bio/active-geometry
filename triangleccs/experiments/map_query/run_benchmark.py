#!/usr/bin/env python3
"""Probe the live Atlas query pipeline through TriangleCCS.

Does not certify Form. Records map behaviour (identify / predict / 129D
tangent) and reads it as a consumer: drop radial when identifiable, ignore
live κ, emit candidate θ against Form anchors.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from triangleccs.datum.form import Form  # noqa: E402
from triangleccs.datum.gauge import wrap_pi  # noqa: E402

from analyze import summarize_run, wrap_atlas_theta  # noqa: E402
from harness import AtlasClient, strip_meta  # noqa: E402
from panel import length_ladder, load_panel, null_sequences  # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def compact_identify(row: dict) -> dict:
    keep = [
        "_ok",
        "_http",
        "lineage",
        "species",
        "genus",
        "family",
        "domain",
        "zone",
        "confidence",
        "prediction_set_size",
        "best_distance",
        "margin",
        "atlas_r",
        "atlas_theta",
        "coordinates",
        "kappa",
        "candidates",
        "encoder_used",
        "latency_ms",
        "error",
    ]
    return {k: row[k] for k in keep if k in row}


def compact_predict(row: dict) -> dict:
    keep = [
        "_ok",
        "_http",
        "domain",
        "domain_confidence",
        "family",
        "family_confidence",
        "genus",
        "genus_confidence",
        "coordinates",
        "radius",
        "kappa",
        "tangent",
        "error",
    ]
    return {k: row[k] for k in keep if k in row}


def pipeline_snapshot(health: dict, info: dict) -> dict:
    return {
        "what_is_wired": {
            "public_globe": "https://www.biosphereatlas.com ball_data.json (v9 consensus, 121351 points)",
            "query_ui": "query-panel.js POST /identify then client-side species lookup in the ball",
            "inference_api": "https://api.biosphereatlas.com FastAPI 15.5.5",
            "encoder": "V15Model dual-path-radial-depth, 129D Poincaré, live GPU",
            "placement": "/identify geodesic hierarchical descent + conformal zones",
            "classify": "/predict domain/family/genus heads + 3D PCA; 129D tangent with API key",
            "batch": "POST /place FASTA",
        },
        "what_is_not_wired": [
            "TriangleCCS Form hash on responses",
            "published 129D→chart transform (examples/atlas_transform.v1.json is an 8D fixture)",
            "sextant (needs a shared aligned panel)",
            "drop-radial on the public path (API still returns atlas_r and live κ)",
            "Address emission (θ candidate / r advisory tags)",
            "freeze-gate",
            "KESTREL short-read path (documented; all observed encoder_used=atlas)",
        ],
        "version_lineage": {
            "docs_map_epoch": "v10.9",
            "live_api": health.get("model_version"),
            "openapi": "15.5.5",
            "ball_index": "v9 consensus",
            "note": "v15.5 is the operational successor of the v10.9 map epoch",
        },
        "live_device": {
            "gpu": health.get("gpu_name"),
            "device": health.get("device"),
            "live_kappa": health.get("kappa"),
            "architecture": info.get("architecture"),
            "prototype_counts": info.get("counts"),
        },
        "constitution": {
            "form_kappa": "CONVENTION — do not overwrite with live κ",
            "radius": "ADVISORY — do not read atlas_r as occupancy",
            "theta": "CANDIDATE until freeze-gate",
            "3d_pca": "decoder view of the map, not the polar chart",
        },
    }


def main() -> int:
    form = Form()
    client = AtlasClient()
    RESULTS.mkdir(exist_ok=True)

    health = client.health()
    client.polite()
    info = client.model_info()
    if not health.get("_ok"):
        print("API not healthy:", health, file=sys.stderr)
        return 1

    panel = load_panel()
    payload: dict = {
        "run_at": datetime.now(timezone.utc).isoformat(),
        "api": client.base,
        "health": strip_meta(health),
        "model_info": strip_meta(info),
        "pipeline": pipeline_snapshot(health, info),
        "form_hash": form.form_hash,
        "panel_meta": [
            {k: v for k, v in p.items() if k != "sequence"} for p in panel
        ],
        "panel": panel,
        "identify": {},
        "predict": {},
    }

    print(f"health {health.get('model_version')} κ_live={health.get('kappa')} gpu={health.get('gpu_name')}")
    print(f"Form {form.version} hash={form.form_hash} κ_convention={form.kappa}")

    for p in panel:
        ident = compact_identify(client.identify(p["sequence"]))
        client.polite()
        pred = compact_predict(client.predict(p["sequence"]))
        client.polite()
        payload["identify"][p["id"]] = ident
        payload["predict"][p["id"]] = pred
        sp = ident.get("species", "?")
        z = ident.get("zone")
        th = ident.get("atlas_theta")
        print(
            f"  {p['id']:16s} true={p['domain']:10s} "
            f"id={ident.get('domain')} {sp} zone={z} "
            f"θ={th} r={ident.get('atlas_r')} "
            f"pred={pred.get('domain')} tan={None if pred.get('tangent') is None else len(pred['tangent'])} "
            f"{ident.get('latency_ms')}ms"
        )

    ecoli = next(p for p in panel if p["id"] == "ecoli_k12")
    payload["length_ladder"] = []
    for rung in length_ladder(ecoli["sequence"]):
        ident = compact_identify(client.identify(rung["sequence"]))
        client.polite()
        payload["length_ladder"].append(
            {"kind": rung["kind"], "length": int(rung["length"]), "identify": ident}
        )
        print(
            f"  ladder {rung['kind']:10s} L={rung['length']:4s} "
            f"{ident.get('species')} zone={ident.get('zone')} "
            f"θ={ident.get('atlas_theta')} r={ident.get('atlas_r')} "
            f"enc={ident.get('encoder_used')}"
        )

    payload["nulls"] = []
    for n in null_sequences(ecoli["sequence"]):
        ident = compact_identify(client.identify(n["sequence"]))
        client.polite()
        payload["nulls"].append({"id": n["id"], "kind": n["kind"], "identify": ident})
        print(
            f"  null {n['id']:16s} {ident.get('domain')} {ident.get('species')} "
            f"zone={ident.get('zone')} conf={ident.get('confidence')} "
            f"d={ident.get('best_distance')}"
        )

    reps = []
    for _ in range(3):
        ident = compact_identify(client.identify(ecoli["sequence"]))
        client.polite()
        reps.append(ident)
    thetas = [wrap_atlas_theta(r["atlas_theta"]) for r in reps if r.get("atlas_theta") is not None]
    coords = [r.get("coordinates") for r in reps]
    payload["reproducibility"] = reps
    payload["reproducibility_summary"] = {
        "n": len(reps),
        "species": [r.get("species") for r in reps],
        "atlas_theta_deg": [round(float(np.degrees(t)), 4) for t in thetas],
        "theta_spread_deg": float(np.degrees(np.ptp(thetas))) if thetas else None,
        "coordinates": coords,
        "same_species": len(set(r.get("species") for r in reps)) == 1,
    }
    print("  repro", payload["reproducibility_summary"])

    analysis = summarize_run(payload, form)
    payload["analysis"] = analysis
    # sequences are large; keep them in panel for analysis then drop from disk dump
    slim = dict(payload)
    slim["panel"] = payload["panel_meta"]
    del slim["panel_meta"]
    if "_transform_full" in slim:
        slim["candidate_transform"] = slim.pop("_transform_full")
    slim.pop("_coords", None)
    # drop tangents from disk? keep them — they are the instrument. ~18*129 numbers is fine.
    out_path = RESULTS / "latest.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(slim, f, indent=2)
    print(f"wrote {out_path}")

    a = analysis
    print("\n=== TriangleCCS reading ===")
    print(f"identify domain accuracy: {a.get('identify_domain_accuracy')}")
    print(f"predict domain accuracy:  {a.get('predict_domain_accuracy')}")
    print(f"zones: {a.get('zones')}")
    print(f"radial head: {a.get('radial_head')}")
    print(f"vector kind: {a.get('vector_kind')} median||v||={a.get('median_vector_norm')}")
    print(f"2D explained: {a.get('candidate_explained_2d')}")
    print(f"atlasθ vs SVD median deg: {a.get('atlas_theta_vs_svd_median_deg')}")
    print(f"Mash vs Poincaré Spearman: {a.get('mash_vs_ambient_poincare_spearman')}")
    print(f"sibling sep deg: {a.get('sibling_median_sep_deg')} distant: {a.get('distant_median_sep_deg')}")
    print(f"freeze-gate: {a.get('freeze_gate')}")
    if a.get("length_ladder"):
        print(f"length ladder θ range deg: {a['length_ladder'].get('theta_range_deg')}")
    print("candidate θ (deg):", a.get("candidate_theta_deg"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
