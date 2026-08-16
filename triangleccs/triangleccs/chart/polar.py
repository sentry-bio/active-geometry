"""Polar chart on H² — θ candidate/certified, r advisory.

Interpreting radius as process depth and angle as divergence is a modelling
choice, not a consequence of the capacity theorem. Origin is the chart origin,
never LUCA.
"""

from __future__ import annotations

import numpy as np

from triangleccs.chart.poincare import logmap0
from triangleccs.datum.form import Form


def euclidean_polar(xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ball coordinates → (r_eucl, θ). r_eucl is advisory chart radius."""
    xy = np.atleast_2d(np.asarray(xy, dtype=np.float64))
    r = np.linalg.norm(xy, axis=1)
    theta = np.arctan2(xy[:, 1], xy[:, 0])
    return r, theta


def from_polar(r: np.ndarray, theta: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=np.float64)
    theta = np.asarray(theta, dtype=np.float64)
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def hyperbolic_radius(xy: np.ndarray, form: Form | None = None) -> np.ndarray:
    """Geodesic radius from the chart origin (still ADVISORY as a depth proxy)."""
    form = form or Form()
    v = logmap0(xy, form.kappa)
    return np.linalg.norm(v, axis=1)
