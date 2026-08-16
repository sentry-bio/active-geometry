"""Register a map realization onto the Form's polar chart.

The map (e.g. Atlas 129D) places genomes; registration emits comparable
Addresses. Radial coordinate remains ADVISORY. The map is a consumer of
Form, not a producer of κ.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from triangleccs.address import Address, make_address
from triangleccs.chart.poincare import distance_matrix, exp_map_zero, logmap0, poincare_distance
from triangleccs.classifier.quartets import classify_quartet_resolvability, measure_quartet_defect
from triangleccs.datum.form import Form
from triangleccs.datum.gauge import wrap_pi
from triangleccs.ledger import Tag
from triangleccs.packing.bound import chart_block_separation


def load_transform(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _backbone(transform: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    if "backbone_basis" in transform:
        basis = np.asarray(transform["backbone_basis"], dtype=np.float64)
        mu = np.asarray(transform["tangent_mean"], dtype=np.float64)
    else:
        basis = np.asarray(transform["v9_backbone_basis_2x129"], dtype=np.float64)
        mu = np.asarray(transform["v9_tangent_mean_129"], dtype=np.float64)
    if basis.shape[0] != 2:
        raise ValueError("backbone_basis must have shape (2, D)")
    return basis, mu


def register_coords(
    atlas_coords: np.ndarray,
    transform: Mapping[str, Any],
    meridian_coord: np.ndarray,
    chirality_coord: np.ndarray | None = None,
    form: Form | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Map coords → (θ, advisory geodesic r on the 2D chart)."""
    th, r, _xy = register_chart(
        atlas_coords, transform, meridian_coord, chirality_coord, form
    )
    return th, r


def register_chart(
    atlas_coords: np.ndarray,
    transform: Mapping[str, Any],
    meridian_coord: np.ndarray,
    chirality_coord: np.ndarray | None = None,
    form: Form | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map coords → (θ, advisory geodesic r, Poincaré xy).

    Radius is the 2D-chart geodesic from the origin (||tangent|| after the
    backbone), tagged ADVISORY. It is not Euclidean ||x|| in the ambient ball.
    """
    form = form or Form()
    coords = np.atleast_2d(np.asarray(atlas_coords, dtype=np.float64))
    basis, mu = _backbone(transform)
    P = (logmap0(coords, form.kappa) - mu) @ basis.T
    th = np.arctan2(P[:, 1], P[:, 0])
    mer = np.asarray(meridian_coord, dtype=np.float64).reshape(1, -1)
    th0 = float(np.arctan2(*((logmap0(mer, form.kappa) - mu) @ basis.T)[0][::-1]))
    th = wrap_pi(th - th0)
    if chirality_coord is not None:
        chir = np.asarray(chirality_coord, dtype=np.float64).reshape(1, -1)
        thm = float(
            np.arctan2(*((logmap0(chir, form.kappa) - mu) @ basis.T)[0][::-1]) - th0
        )
        if wrap_pi(np.array([thm]))[0] < 0:
            th = -th
    r = np.linalg.norm(P, axis=1)
    xy = exp_map_zero(np.column_stack([r * np.cos(th), r * np.sin(th)]), form.kappa)
    return th, r, xy


def _projection_residual(
    atlas_coords: np.ndarray, xy: np.ndarray, form: Form
) -> np.ndarray:
    """Per-point RMSE of 2D chart distances vs ambient Poincaré distances."""
    coords = np.atleast_2d(np.asarray(atlas_coords, dtype=np.float64))
    n = coords.shape[0]
    chart = distance_matrix(xy, form)
    out = np.zeros(n, dtype=np.float64)
    if n < 2:
        return out
    ambient = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        ambient[i] = poincare_distance(coords[i], coords, form.kappa)
    np.fill_diagonal(ambient, 0.0)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        err = chart[i, mask] - ambient[i, mask]
        out[i] = float(np.sqrt(np.mean(err * err)))
    return out


def addresses_from_registration(
    atlas_coords: np.ndarray,
    transform: Mapping[str, Any],
    meridian_coord: np.ndarray,
    chirality_coord: np.ndarray | None = None,
    form: Form | None = None,
    *,
    theta_status: str = "candidate",
    max_quartets: int = 2_000,
    seed: int = 0,
) -> list[Address]:
    form = form or Form()
    th, r, xy = register_chart(
        atlas_coords, transform, meridian_coord, chirality_coord, form
    )
    extra = {}
    if transform.get("inheritance") == "warm_start":
        extra["inheritance"] = Tag.CIRCULAR.value
    block = float(chart_block_separation(xy, form))
    residual = _projection_residual(atlas_coords, xy, form)
    n = len(th)
    if n >= 4:
        chart_d = distance_matrix(xy, form)
        defect = measure_quartet_defect(chart_d, max_quartets=max_quartets, seed=seed)
        reso = classify_quartet_resolvability(
            chart_d, jc_ceiling=None, max_quartets=max_quartets, seed=seed
        )
        delta = float(defect["delta_q"]["q50"])
        resolvable = float(reso["resolvable_fraction"])
    else:
        delta = 0.0
        resolvable = 1.0
    out: list[Address] = []
    for i in range(n):
        out.append(
            make_address(
                form=form,
                theta=float(th[i]),
                r=float(r[i]),
                theta_status=theta_status,  # type: ignore[arg-type]
                delta=delta,
                resolvable=resolvable,
                block_sep=block,
                residual=float(residual[i]),
                extra_tags=extra or None,
            )
        )
    return out


def cross_instrument_agreement(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    meridian_index: int,
    chirality_index: int,
    form: Form | None = None,
) -> dict[str, float | int | str]:
    """Angular median residual between two instruments (θ chart agreement)."""
    from triangleccs.chart.polar import hyperbolic_radius
    from triangleccs.datum.gauge import svd_backbone_theta

    form = form or Form()
    th_a = svd_backbone_theta(coords_a, meridian_index, chirality_index, form)
    th_b = svd_backbone_theta(coords_b, meridian_index, chirality_index, form)
    d = np.abs(wrap_pi(th_a - th_b))
    ra = hyperbolic_radius(np.asarray(coords_a, dtype=np.float64), form)
    rb = hyperbolic_radius(np.asarray(coords_b, dtype=np.float64), form)
    order = lambda x: np.argsort(np.argsort(x))
    rho = float(np.corrcoef(order(ra), order(rb))[0, 1])
    return {
        "angular_median_deg": float(np.degrees(np.median(d))),
        "within_30deg_frac": float((d < np.pi / 6).mean()),
        "radial_advisory_spearman": rho,
        "n": int(len(th_a)),
        "note": (
            "Chart agreement within the imposed model; radial is ADVISORY geodesic; "
            "not evidence of information-tightness or saturation."
        ),
    }
