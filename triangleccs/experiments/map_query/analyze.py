"""TriangleCCS reading of map outputs.

Form κ is CONVENTION (1.25). Live API κ is a map field and is not copied
onto Form. Radius is ADVISORY. θ from SVD of the 129D tangent (radial head
dropped when identifiable) is CANDIDATE. Mash k-mer distance is a second
map, not a sextant: the sequences are unaligned.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from triangleccs.chart.poincare import ball_radius, exp_map_zero, logmap0, poincare_distance
from triangleccs.datum.form import Form
from triangleccs.datum.gauge import wrap_pi
from triangleccs.datum.registration import addresses_from_registration

from panel import canonical_kmers, mash_distance


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    if np.std(ra) == 0 or np.std(rb) == 0:
        return float("nan")
    return float(np.corrcoef(ra.astype(float), rb.astype(float))[0, 1])


def interpret_vectors(vectors: np.ndarray, form: Form) -> dict[str, Any]:
    V = np.atleast_2d(np.asarray(vectors, dtype=np.float64))
    norms = np.linalg.norm(V, axis=1)
    br = ball_radius(form.kappa)
    median = float(np.median(norms))
    if median < br * 0.98:
        kind = "poincare_coords"
        coords = V
    else:
        kind = "tangent"
        coords = exp_map_zero(V, form.kappa)
    return {
        "kind": kind,
        "median_norm": median,
        "ball_radius_form": br,
        "coords": coords,
        "norms": norms.tolist(),
    }


def drop_radial_head(
    vectors: np.ndarray, radius_obs: np.ndarray
) -> dict[str, Any]:
    V = np.atleast_2d(np.asarray(vectors, dtype=np.float64))
    r = np.asarray(radius_obs, dtype=np.float64)
    corrs = []
    for j in range(V.shape[1]):
        col = V[:, j]
        if np.std(col) < 1e-12 or np.std(r) < 1e-12:
            corrs.append(0.0)
        else:
            corrs.append(float(np.corrcoef(col, r)[0, 1]))
    abs_c = np.abs(np.asarray(corrs))
    j = int(np.argmax(abs_c))
    dropped = abs_c[j] >= 0.5
    kept = np.delete(V, j, axis=1) if dropped else V
    return {
        "n_dim": int(V.shape[1]),
        "radial_index": j,
        "radial_corr": float(corrs[j]),
        "dropped": bool(dropped),
        "kept_dim": int(kept.shape[1]),
        "all_abs_corr_max": float(abs_c.max()),
        "last_dim_corr": float(corrs[-1]) if corrs else None,
        "vectors": kept,
        "corrs": [float(c) for c in corrs],
    }


def candidate_backbone(
    coords: np.ndarray,
    meridian_index: int,
    chirality_index: int,
    form: Form,
) -> dict[str, Any]:
    T = logmap0(coords, form.kappa)
    mu = T.mean(axis=0)
    _, svals, vt = np.linalg.svd(T - mu, full_matrices=False)
    p = (T - mu) @ vt[:2].T
    th = wrap_pi(np.arctan2(p[:, 1], p[:, 0]) - np.arctan2(p[meridian_index, 1], p[meridian_index, 0]))
    if th[chirality_index] < 0:
        th = -th
        vt = vt.copy()
        vt[1] *= -1
        p = (T - mu) @ vt[:2].T
        th = wrap_pi(np.arctan2(p[:, 1], p[:, 0]) - np.arctan2(p[meridian_index, 1], p[meridian_index, 0]))
    r = np.linalg.norm(p, axis=1)
    explained = (svals[:2] ** 2).sum() / max((svals ** 2).sum(), 1e-12)
    transform = {
        "form_version": form.version,
        "form_hash": form.form_hash,
        "map_id": "v15.5-live-16s-svd-candidate",
        "backbone_basis": vt[:2].tolist(),
        "tangent_mean": mu.tolist(),
        "certified": False,
        "inheritance": "independent",
        "note": (
            "Candidate SVD backbone on this run's marker panel. "
            "Not a published transform. Freeze-gate cannot pass without a "
            "sextant-placed aligned panel."
        ),
    }
    addrs = addresses_from_registration(
        coords, transform, coords[meridian_index], coords[chirality_index], form=form
    )
    return {
        "theta": [float(x) for x in th],
        "r_advisory": [float(x) for x in r],
        "explained_2d": float(explained),
        "singular_values_head": [float(x) for x in svals[:8]],
        "transform": transform,
        "addresses": [a.to_dict() for a in addrs],
    }


def wrap_atlas_theta(theta: float) -> float:
    return float(wrap_pi(np.array([theta]))[0])


def domain_angular_stats(theta: np.ndarray, domains: list[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for d in sorted(set(domains)):
        sub = theta[np.array(domains) == d]
        out[d] = {
            "n": int(len(sub)),
            "median_deg": float(np.degrees(np.median(sub))) if len(sub) else None,
            "span_deg": float(np.degrees(np.ptp(sub))) if len(sub) else None,
        }
    # pairwise circular median separation of domain medians
    med = {d: float(np.median(theta[np.array(domains) == d])) for d in set(domains)}
    pairs = {}
    names = sorted(med)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            pairs[f"{a}|{b}"] = float(np.degrees(abs(wrap_pi(np.array([med[a] - med[b]]))[0])))
    out["median_separations_deg"] = pairs
    return out


def mash_matrix(sequences: list[str], k: int = 21) -> np.ndarray:
    sets = [canonical_kmers(s, k=k) for s in sequences]
    n = len(sets)
    d = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d[i, j] = d[j, i] = mash_distance(sets[i], sets[j], k=k)
    return d


def ambient_poincare_matrix(coords: np.ndarray, form: Form) -> np.ndarray:
    c = np.asarray(coords, dtype=np.float64)
    n = c.shape[0]
    out = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        out[i] = poincare_distance(c[i], c, form.kappa)
    np.fill_diagonal(out, 0.0)
    return out


def summarize_run(payload: dict[str, Any], form: Form) -> dict[str, Any]:
    panel = payload["panel"]
    ids = [p["id"] for p in panel]
    domains = [p["domain"] for p in panel]
    seqs = [p["sequence"] for p in panel]
    ident = payload["identify"]
    pred = payload["predict"]

    ident_ok = [ident[i] for i in ids]
    pred_ok = [pred[i] for i in ids]
    domain_id_match = [
        (ident[i].get("domain") == panel[k]["domain"])
        for k, i in enumerate(ids)
        if ident[i].get("_ok")
    ]
    domain_pred_match = [
        (pred[i].get("domain") == panel[k]["domain"])
        for k, i in enumerate(ids)
        if pred[i].get("_ok")
    ]

    atlas_theta = np.array(
        [wrap_atlas_theta(ident[i]["atlas_theta"]) for i in ids if ident[i].get("atlas_theta") is not None],
        dtype=np.float64,
    )
    atlas_r = np.array([ident[i].get("atlas_r") for i in ids], dtype=np.float64)
    zones = [ident[i].get("zone") for i in ids]
    species = [ident[i].get("species") for i in ids]
    encoder = [ident[i].get("encoder_used") for i in ids]

    tangents = []
    tangent_ids = []
    for i, p in zip(ids, pred_ok):
        t = p.get("tangent")
        if t:
            tangents.append(t)
            tangent_ids.append(i)
    analysis: dict[str, Any] = {
        "form": form.summary(),
        "live_kappa": payload["health"].get("kappa"),
        "kappa_note": "live κ is a map field; Form.kappa is CONVENTION and is not updated",
        "n_panel": len(panel),
        "identify_domain_accuracy": float(np.mean(domain_id_match)) if domain_id_match else None,
        "predict_domain_accuracy": float(np.mean(domain_pred_match)) if domain_pred_match else None,
        "encoder_used": {e: encoder.count(e) for e in sorted(set(encoder))},
        "zones": {z: zones.count(z) for z in sorted(set(z for z in zones if z))},
        "species_calls": {i: ident[i].get("species") for i in ids},
        "atlas_theta_domain": domain_angular_stats(
            np.array([wrap_atlas_theta(ident[i]["atlas_theta"]) for i in ids], dtype=np.float64),
            domains,
        ),
        "mean_atlas_r": float(np.nanmean(atlas_r)),
        "reproducibility": payload.get("reproducibility_summary"),
    }

    if len(tangents) >= 4:
        T = np.asarray(tangents, dtype=np.float64)
        interp = interpret_vectors(T, form)
        coords = interp["coords"]
        radial = drop_radial_head(T, atlas_r)
        analysis["vector_kind"] = interp["kind"]
        analysis["median_vector_norm"] = interp["median_norm"]
        analysis["radial_head"] = {
            "index": radial["radial_index"],
            "corr": radial["radial_corr"],
            "dropped": radial["dropped"],
            "kept_dim": radial["kept_dim"],
            "last_dim_corr": radial["last_dim_corr"],
        }
        mer = ids.index("ecoli_k12")
        chir = ids.index("mjannaschii")
        # Prefer angular-only coords when radial head was dropped from tangent
        if radial["dropped"] and interp["kind"] == "tangent":
            ang_coords = exp_map_zero(radial["vectors"], form.kappa)
        else:
            ang_coords = coords
        backbone = candidate_backbone(ang_coords, mer, chir, form)
        analysis["candidate_explained_2d"] = backbone["explained_2d"]
        analysis["candidate_theta_deg"] = {
            ids[i]: round(float(np.degrees(backbone["theta"][i])), 2) for i in range(len(ids))
        }
        analysis["candidate_theta_domain"] = domain_angular_stats(
            np.asarray(backbone["theta"]), domains
        )
        analysis["candidate_transform_certified"] = False
        mash = mash_matrix(seqs)
        ambient = ambient_poincare_matrix(ang_coords, form)
        iu = np.triu_indices(len(ids), k=1)
        analysis["mash_vs_ambient_poincare_spearman"] = _spearman(mash[iu], ambient[iu])
        # atlas_theta vs SVD θ (circular residual)
        at = np.array([wrap_atlas_theta(ident[i]["atlas_theta"]) for i in ids])
        dth = np.abs(wrap_pi(at - np.asarray(backbone["theta"])))
        analysis["atlas_theta_vs_svd_median_deg"] = float(np.degrees(np.median(dth)))
        analysis["atlas_theta_vs_svd_within_30deg"] = float((dth < np.pi / 6).mean())
        analysis["addresses"] = backbone["addresses"]
        analysis["transform"] = {
            k: v
            for k, v in backbone["transform"].items()
            if k not in {"backbone_basis", "tangent_mean"}
        }
        analysis["transform"]["backbone_basis_shape"] = [
            2,
            len(backbone["transform"]["tangent_mean"]),
        ]
        payload["_transform_full"] = backbone["transform"]
        payload["_coords"] = ang_coords.tolist()

        enterobac = [i for i, p in enumerate(panel) if p["clade"] == "Enterobacteriaceae"]
        distant = [i for i, p in enumerate(panel) if p["role"] == "distant_bacteria"]
        if len(enterobac) >= 2 and distant:
            th = np.asarray(backbone["theta"])
            sib = []
            for a in range(len(enterobac)):
                for b in range(a + 1, len(enterobac)):
                    sib.append(abs(wrap_pi(np.array([th[enterobac[a]] - th[enterobac[b]]]))[0]))
            far = []
            e0 = enterobac[0]
            for j in distant:
                far.append(abs(wrap_pi(np.array([th[e0] - th[j]]))[0]))
            analysis["sibling_median_sep_deg"] = float(np.degrees(np.median(sib)))
            analysis["distant_median_sep_deg"] = float(np.degrees(np.median(far)))

        analysis["freeze_gate"] = {
            "evaluable": False,
            "reason": (
                "No aligned sextant panel; SVD θ vs map atlas_theta is two "
                "readings of the same map, not an independent witness. "
                "theta_status remains candidate."
            ),
            "theta_status": "candidate",
            "passed": False,
        }

    # Length ladder: θ stability
    ladder = payload.get("length_ladder") or []
    if ladder:
        thetas = [wrap_atlas_theta(row["identify"]["atlas_theta"]) for row in ladder if row["identify"].get("atlas_theta") is not None]
        rs = [row["identify"].get("atlas_r") for row in ladder]
        analysis["length_ladder"] = {
            "n": len(ladder),
            "atlas_theta_deg": [round(float(np.degrees(t)), 2) for t in thetas],
            "atlas_r": rs,
            "theta_range_deg": float(np.degrees(np.ptp(thetas))) if thetas else None,
            "r_range": float(np.nanmax(rs) - np.nanmin(rs)) if rs else None,
            "species": [row["identify"].get("species") for row in ladder],
            "zone": [row["identify"].get("zone") for row in ladder],
            "encoder": [row["identify"].get("encoder_used") for row in ladder],
            "note": "θ should be more stable than r if the map has a usable backbone; r is ADVISORY",
        }

    nulls = payload.get("nulls") or []
    if nulls:
        analysis["nulls"] = [
            {
                "id": n["id"],
                "zone": n["identify"].get("zone"),
                "species": n["identify"].get("species"),
                "confidence": n["identify"].get("confidence"),
                "best_distance": n["identify"].get("best_distance"),
                "domain": n["identify"].get("domain"),
            }
            for n in nulls
        ]

    analysis["what_this_does_not_claim"] = [
        "Form certification (freeze-gate not evaluable)",
        "live κ as a measurement of curvature",
        "atlas_r as occupancy of exponential room",
        "taxonomy accuracy as occupancy or saturation",
        "3D PCA globe as the polar chart",
    ]
    return analysis
