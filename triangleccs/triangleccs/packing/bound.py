"""Metric packing bound — Layer I limit, curvature-free.

Finite-radius packing counts and occupancy slopes. Occupancy is never promoted
to asymptotic h_pack or to efficiency η. No path from (h, κ) to n.

The bound is stated for an arbitrary pointed metric. Euclidean coordinates of
a Poincaré chart are the wrong metric for that chart; use PoincareMetric.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from triangleccs.metric import EuclideanMetric, Metric, PoincareMetric, metric_from_form


def packing_count(
    points: np.ndarray,
    origin: np.ndarray,
    radius: float,
    epsilon: float,
    metric: Metric | None = None,
) -> int:
    """Cardinality of an ε-separated subset of points inside B(origin, radius).

    Greedy packing of the given finite sample (lower bound on the true packing
    number of the continuum ball). Distances come from ``metric``.
    """
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    if radius < 0:
        raise ValueError("radius must be non-negative")
    metric = metric or EuclideanMetric()
    pts = np.asarray(points, dtype=np.float64)
    o = np.asarray(origin, dtype=np.float64).reshape(-1)
    if pts.ndim != 2:
        raise ValueError("points must have shape (N, d)")
    dists = metric.distance_to(o, pts)
    inside = pts[dists <= radius + 1e-12]
    if inside.size == 0:
        return 0
    dmat = metric.pairwise(inside)
    selected: list[int] = []
    remaining = list(range(inside.shape[0]))
    while remaining:
        i = remaining.pop(0)
        selected.append(i)
        remaining = [
            j for j in remaining if float(dmat[i, j]) >= epsilon - 1e-12
        ]
    return len(selected)


def packing_monotone(
    points: np.ndarray,
    origin: np.ndarray,
    radii: Sequence[float],
    epsilon: float,
    metric: Metric | None = None,
) -> list[int]:
    """Packing counts at increasing radii — must be nondecreasing."""
    metric = metric or EuclideanMetric()
    return [
        packing_count(points, origin, float(r), epsilon, metric=metric)
        for r in radii
    ]


def chart_packing_count(
    points: np.ndarray,
    radius: float,
    epsilon: float,
    form: object | None = None,
) -> int:
    """Packing in the Form's Poincaré metric, origin at the chart origin."""
    metric = metric_from_form(form)
    dim = int(np.asarray(points, dtype=np.float64).shape[1])
    origin = np.zeros(dim, dtype=np.float64)
    return packing_count(points, origin, radius, epsilon, metric=metric)


def fit_ball_occupancy_slope(
    distances: Sequence[float],
    *,
    min_count: int = 3,
    max_fraction: float = 0.5,
    points: int = 32,
) -> dict[str, float | int | str] | None:
    """Finite-sample log-occupancy vs radius slope. Not asymptotic h_pack.

    ``distances`` are already in the caller's metric (Euclidean or Poincaré).
    """
    ordered = np.sort(np.asarray(distances, dtype=np.float64))
    n = ordered.size
    max_count = min(n - 1, int(math.floor(max_fraction * n)))
    if max_count - min_count < 2:
        return None
    requested = np.unique(
        np.linspace(min_count, max_count, num=min(points, max_count), dtype=int)
    )
    radii = ordered[np.clip(requested - 1, 0, n - 1)]
    occupancies = np.searchsorted(ordered, radii, side="right")
    valid = radii > 0
    radii = radii[valid]
    occupancies = occupancies[valid]
    unique_radii, indices = np.unique(radii, return_index=True)
    occupancies = occupancies[indices]
    if unique_radii.size < 3:
        return None
    radius_scale = float(unique_radii.max())
    scaled = unique_radii / radius_scale
    design = np.column_stack([np.ones(unique_radii.size), scaled])
    response = np.log(occupancies.astype(np.float64))
    intercept, scaled_slope = np.linalg.lstsq(design, response, rcond=None)[0]
    prediction = intercept + scaled_slope * scaled
    slope = float(scaled_slope / radius_scale)
    residual = float(np.sum((response - prediction) ** 2))
    total = float(np.sum((response - response.mean()) ** 2))
    r_squared = 1.0 - residual / total if total > 0 else 1.0
    return {
        "finite_ball_occupancy_slope": max(0.0, slope),
        "fit_r_squared": r_squared,
        "fit_points": int(unique_radii.size),
        "status": "INSTRUMENT",
        "note": "finite-sample occupancy slope; not h_pack; not η",
    }


def block_separation_fraction(
    points: np.ndarray,
    epsilon: float,
    metric: Metric | None = None,
) -> float:
    """Fraction of points that remain after greedy ε-separation of the whole set."""
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[0] == 0:
        return 0.0
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    metric = metric or EuclideanMetric()
    dmat = metric.pairwise(pts)
    n = pts.shape[0]
    selected: list[int] = []
    for i in range(n):
        if all(float(dmat[i, j]) >= epsilon - 1e-12 for j in selected):
            selected.append(i)
    return len(selected) / n


def chart_block_separation(
    points: np.ndarray,
    form: object | None = None,
) -> float:
    """Block-separation fraction in the Form's Poincaré metric at Form.epsilon."""
    from triangleccs.datum.form import Form

    f = form if isinstance(form, Form) else Form()
    return block_separation_fraction(
        points, f.epsilon, metric=PoincareMetric(f.kappa)
    )
