"""Two-anchor O(2) gauge: meridian + chirality.

Prime meridian (E. coli) fixes θ = 0. Chirality anchor (M. jannaschii) fixes
reflection. Leftover uniqueness is global O(2) until both anchors are applied.
"""

from __future__ import annotations

import numpy as np

from triangleccs.chart.poincare import logmap0
from triangleccs.datum.form import Form


def wrap_pi(theta: np.ndarray) -> np.ndarray:
    return (theta + np.pi) % (2.0 * np.pi) - np.pi


def anchor_theta(
    xy: np.ndarray,
    meridian_xy: np.ndarray,
    chirality_xy: np.ndarray | None = None,
    form: Form | None = None,
) -> np.ndarray:
    """Return θ in radians with meridian at 0 and optional chirality fix."""
    form = form or Form()
    xy = np.atleast_2d(np.asarray(xy, dtype=np.float64))
    # Work in tangent space so the angular chart is geometrically honest.
    T = logmap0(xy, form.kappa)
    t0 = logmap0(np.asarray(meridian_xy, dtype=np.float64).reshape(1, -1), form.kappa)[0]
    th = np.arctan2(T[:, 1], T[:, 0])
    th0 = float(np.arctan2(t0[1], t0[0]))
    th = wrap_pi(th - th0)
    if chirality_xy is not None:
        tc = logmap0(
            np.asarray(chirality_xy, dtype=np.float64).reshape(1, -1), form.kappa
        )[0]
        thc = wrap_pi(float(np.arctan2(tc[1], tc[0])) - th0)
        if thc < 0:
            th = -th
    return th


def svd_backbone_theta(
    coords: np.ndarray,
    meridian_index: int,
    chirality_index: int,
    form: Form | None = None,
) -> np.ndarray:
    """Gauge-free-ish 2D backbone θ from any ambient Poincaré coords (N, D)."""
    form = form or Form()
    T = logmap0(coords, form.kappa)
    mu = T.mean(axis=0)
    _, _, vt = np.linalg.svd(T - mu, full_matrices=False)
    p = (T - mu) @ vt[:2].T
    th = np.arctan2(p[:, 1], p[:, 0])
    th = wrap_pi(th - th[meridian_index])
    if th[chirality_index] < 0:
        th = -th
    return th
