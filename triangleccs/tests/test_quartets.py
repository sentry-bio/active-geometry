"""Quartet classifier."""

from __future__ import annotations

import numpy as np

from triangleccs.classifier.quartets import (
    measure_quartet_defect,
    path_tree_distance,
    quartet_split,
)


def test_exact_tree_delta_zero():
    d = path_tree_distance(8, branch=1.0)
    out = measure_quartet_defect(d, max_quartets=70, seed=0)
    assert out["exact_fraction"] == 1.0
    assert out["delta_q"]["q95"] == 0.0


def test_noise_not_exact():
    rng = np.random.default_rng(1)
    n = 10
    a = rng.random((n, n))
    d = (a + a.T) / 2.0
    np.fill_diagonal(d, 0.0)
    # enforce positivity
    d = d + 0.1
    np.fill_diagonal(d, 0.0)
    out = measure_quartet_defect(d, max_quartets=200, seed=0)
    assert out["exact_fraction"] < 0.5


def test_quartet_split_star_is_tied_but_defined():
    d = path_tree_distance(4, branch=1.0)
    # Star: all three pair-sums equal 4; min() is stable on insertion order.
    split = quartet_split(d, 0, 1, 2, 3)
    assert split == frozenset({frozenset({0, 1}), frozenset({2, 3})})
