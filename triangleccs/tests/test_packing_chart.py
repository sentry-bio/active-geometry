"""Packing and chart geometry."""

from __future__ import annotations

import numpy as np

from triangleccs.chart.poincare import distance_matrix, exp_map_zero, logmap0, poincare_distance
from triangleccs.datum.form import Form
from triangleccs.datum.gauge import anchor_theta
from triangleccs.metric import EuclideanMetric, PoincareMetric
from triangleccs.packing.bound import (
    chart_packing_count,
    packing_count,
    packing_monotone,
)


def test_packing_monotone():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(40, 2)) * 0.2
    origin = np.zeros(2)
    radii = [0.1, 0.2, 0.5, 1.0]
    counts = packing_monotone(pts, origin, radii, epsilon=0.05)
    assert counts == sorted(counts)


def test_packing_count_positive():
    pts = np.array([[0.0, 0.0], [0.2, 0.0], [0.0, 0.2]], dtype=float)
    assert packing_count(pts, np.zeros(2), 1.0, 0.1) >= 1


def test_poincare_packing_disagrees_with_euclidean_near_boundary():
    """Limit and chart must share a metric: Euclidean ||x-y|| is the wrong one."""
    form = Form()
    br = form.ball_radius
    # Two points near the ball boundary, close in Euclidean coordinates.
    x = 0.92 * br
    pts = np.array([[x, 0.0], [x, 0.006], [0.0, 0.0]], dtype=float)
    origin = np.zeros(2)
    eu = packing_count(
        pts, origin, radius=5.0, epsilon=0.05, metric=EuclideanMetric()
    )
    po = packing_count(
        pts, origin, radius=5.0, epsilon=0.05, metric=PoincareMetric(form.kappa)
    )
    d_eu = float(np.linalg.norm(pts[0] - pts[1]))
    d_po = float(poincare_distance(pts[0], pts[1], form.kappa))
    assert d_eu < 0.05 < d_po
    assert eu < po
    assert chart_packing_count(pts, radius=5.0, epsilon=0.05, form=form) == po


def test_poincare_packing_monotone():
    form = Form()
    rng = np.random.default_rng(1)
    v = rng.normal(size=(24, 2)) * 0.15
    pts = exp_map_zero(v, form.kappa)
    counts = packing_monotone(
        pts,
        np.zeros(2),
        [0.05, 0.2, 0.6, 1.5],
        epsilon=0.05,
        metric=PoincareMetric(form.kappa),
    )
    assert counts == sorted(counts)


def test_poincare_roundtrip_origin():
    form = Form()
    v = np.array([[0.3, -0.1], [0.0, 0.2]])
    x = exp_map_zero(v, form.kappa)
    back = logmap0(x, form.kappa)
    assert np.allclose(back, v, atol=1e-5)


def test_gauge_meridian_zero():
    form = Form()
    pts = exp_map_zero(np.array([[0.4, 0.0], [0.0, 0.3], [-0.2, 0.1]]), form.kappa)
    th = anchor_theta(pts, pts[0], pts[1], form)
    assert abs(th[0]) < 1e-9


def test_distance_matrix_symmetric():
    form = Form()
    pts = exp_map_zero(np.array([[0.2, 0.0], [0.0, 0.2], [0.1, 0.1]]), form.kappa)
    d = distance_matrix(pts, form)
    assert np.allclose(d, d.T)
    assert np.allclose(np.diag(d), 0.0)
