"""Poincaré ball model of H² at a frozen κ.

Theorem 4.4 of Active Geometry (paper; Skenderi lower bound) establishes
existence: weighted relational capacity of real hyperbolic space equals block
capacity at exponential order. That is why this host may be inhabited without
an exponential genealogical tax. It is not a uniqueness theorem (curvature
genericity remains an overlay) and it does not license κ = (h ln 2)².
"""

from __future__ import annotations

import math

import numpy as np

from triangleccs.datum.form import Form


def ball_radius(kappa: float) -> float:
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    return 1.0 / math.sqrt(kappa)


def project_to_ball(x: np.ndarray, kappa: float, eps: float = 1e-5) -> np.ndarray:
    """Project onto the open ball {x : κ||x||² < 1}."""
    x = np.asarray(x, dtype=np.float64)
    max_norm = ball_radius(kappa) - eps
    norms = np.linalg.norm(x, axis=-1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    scale = np.minimum(1.0, max_norm / norms)
    return x * scale


def exp_map_zero(v: np.ndarray, kappa: float) -> np.ndarray:
    """Exponential map at the origin of the Poincaré ball."""
    v = np.asarray(v, dtype=np.float64)
    sqrt_c = math.sqrt(kappa)
    v_norm = np.linalg.norm(v, axis=-1, keepdims=True)
    v_norm = np.maximum(v_norm, 1e-12)
    return np.tanh(sqrt_c * v_norm / 2.0) * v / (sqrt_c * v_norm)


def logmap0(x: np.ndarray, kappa: float) -> np.ndarray:
    """Poincaré ball → tangent space at the origin."""
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    br = ball_radius(kappa)
    n = np.linalg.norm(x, axis=1, keepdims=True)
    n = np.clip(n, 1e-12, br * (1.0 - 1e-6))
    return (2.0 / math.sqrt(kappa)) * np.arctanh(math.sqrt(kappa) * n) * (x / n)


def poincare_distance(
    x: np.ndarray, y: np.ndarray, kappa: float
) -> np.ndarray:
    """Geodesic distance in the Poincaré ball with curvature −κ."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    sqrt_c = math.sqrt(kappa)
    x_sq = np.clip((x * x).sum(axis=-1), 0.0, 1.0 / kappa - 1e-5)
    y_sq = np.clip((y * y).sum(axis=-1), 0.0, 1.0 / kappa - 1e-5)
    diff_sq = ((x - y) ** 2).sum(axis=-1)
    num = 2.0 * kappa * diff_sq
    denom = (1.0 - kappa * x_sq) * (1.0 - kappa * y_sq)
    arg = 1.0 + num / np.maximum(denom, 1e-12)
    return (1.0 / sqrt_c) * np.arccosh(np.maximum(arg, 1.0 + 1e-7))


def distance_matrix(points: np.ndarray, form: Form | None = None) -> np.ndarray:
    """Pairwise Poincaré distances for an (N, 2) chart cloud."""
    form = form or Form()
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("chart points must have shape (N, 2)")
    n = pts.shape[0]
    out = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        out[i] = poincare_distance(pts[i], pts, form.kappa)
    np.fill_diagonal(out, 0.0)
    return out
