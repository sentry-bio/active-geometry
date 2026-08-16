"""Sanity checks for the E5 small-scale runner (generators and geometry)."""

from __future__ import annotations

import math
import unittest

import numpy as np

from experiments.e5_trained_hierarchy import (
    complete_tree,
    euclidean_pairwise,
    grid_distances,
    poincare_pairwise,
    poincare_radius,
    to_ball,
)


class E5GeometryTests(unittest.TestCase):
    def test_binary_tree_root_degree_and_leaf_count(self) -> None:
        dist, depths = complete_tree(2, 3)
        self.assertEqual(dist.shape, (15, 15))
        self.assertEqual(int(np.sum(depths == 3)), 8)
        self.assertEqual(dist[0, 1], 1.0)
        self.assertEqual(dist[1, 2], 2.0)

    def test_grid_is_polynomial_l1(self) -> None:
        dist, depths = grid_distances(4)
        self.assertEqual(dist.shape, (16, 16))
        self.assertEqual(dist[0, -1], 6.0)
        self.assertEqual(depths[0], 0.0)

    def test_poincare_origin_radius_matches_formula(self) -> None:
        origin = np.zeros((1, 2))
        p = np.array([[0.5, 0.0]])
        d = poincare_pairwise(np.vstack([origin, p]))[0, 1]
        self.assertAlmostEqual(d, poincare_radius(p)[0], places=6)

    def test_to_ball_stays_inside(self) -> None:
        z = np.array([[10.0, 10.0], [0.0, 0.0]])
        p = to_ball(z)
        self.assertLess(np.linalg.norm(p[0]), 1.0)
        self.assertEqual(np.linalg.norm(p[1]), 0.0)

    def test_euclidean_pairwise_agrees_with_norm(self) -> None:
        pts = np.array([[0.0, 0.0], [3.0, 4.0]])
        self.assertAlmostEqual(euclidean_pairwise(pts)[0, 1], 5.0)


if __name__ == "__main__":
    unittest.main()
