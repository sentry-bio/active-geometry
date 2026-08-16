"""Pointed metrics — packing is curvature-free; the chart supplies one metric.

Layer I counts ε-separated points in a pointed metric space. Euclidean
coordinates of a Poincaré chart are not that metric. Callers that inhabit
the Form must pass ``PoincareMetric(form.kappa)``.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from triangleccs.chart.poincare import poincare_distance


class Metric(Protocol):
    def distance_to(self, origin: np.ndarray, points: np.ndarray) -> np.ndarray:
        """Geodesic radius of each point from origin, shape (N,)."""

    def pairwise(self, points: np.ndarray) -> np.ndarray:
        """Pairwise distances, shape (N, N)."""


class EuclideanMetric:
    """Flat metric on ambient coordinates. Used for Layer I unit tests."""

    def distance_to(self, origin: np.ndarray, points: np.ndarray) -> np.ndarray:
        o = np.asarray(origin, dtype=np.float64).reshape(-1)
        pts = np.asarray(points, dtype=np.float64)
        return np.linalg.norm(pts - o, axis=1)

    def pairwise(self, points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        diff = pts[:, None, :] - pts[None, :, :]
        return np.linalg.norm(diff, axis=-1)


class PoincareMetric:
    """Geodesic metric of the Poincaré ball at frozen κ."""

    def __init__(self, kappa: float) -> None:
        if kappa <= 0:
            raise ValueError("kappa must be positive")
        self.kappa = float(kappa)

    def distance_to(self, origin: np.ndarray, points: np.ndarray) -> np.ndarray:
        o = np.asarray(origin, dtype=np.float64).reshape(-1)
        pts = np.asarray(points, dtype=np.float64)
        return np.asarray(poincare_distance(o, pts, self.kappa), dtype=np.float64)

    def pairwise(self, points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        n = pts.shape[0]
        out = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            out[i] = poincare_distance(pts[i], pts, self.kappa)
        np.fill_diagonal(out, 0.0)
        return out


def metric_from_form(form: object | None = None) -> Metric:
    """Chart metric. ``form`` is ``Form``; default Form() if omitted."""
    from triangleccs.datum.form import Form

    f = form if isinstance(form, Form) else Form()
    return PoincareMetric(f.kappa)
