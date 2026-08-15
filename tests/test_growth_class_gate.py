"""Layer 0: identities, measurability predicate, and the 2×-span death."""

from __future__ import annotations

import math
import unittest

import numpy as np

from tests.test_addressability_meter import balanced_binary_tree_distances
from tools.growth_class_gate import (
    Z_DETECT,
    analyze_distances,
    classify_occupancy,
    critical_sample_size,
    hellinger_affinity,
    lecam_risk_lower_bound,
    measurability,
    midpoint_exponent,
    midpoint_ratio,
    occupancy_profile,
    restrict_span,
    small_span_leading,
    span_information,
)


def integer_grid_distances(side: int) -> np.ndarray:
    """ℓ¹ distances on a 2-D integer grid — polynomial host."""
    coords = np.array([(i, j) for i in range(side) for j in range(side)], dtype=float)
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sum(np.abs(diff), axis=-1)


def log_gap(d: float, r: float, t: float) -> float:
    return d * math.log(t) - d * (math.log(r) / (r - 1.0)) * (t - 1.0)


class IdentityTests(unittest.TestCase):
    def test_midpoint_identity(self) -> None:
        for r in (1.5, 2.0, math.e, 3.0, 4.0, 8.0):
            for d in (1.0, 2.0, 3.0):
                self.assertAlmostEqual(
                    log_gap(d, r, math.sqrt(r)),
                    midpoint_exponent(r, d),
                    places=12,
                )

    def test_max_gap_at_log_mean(self) -> None:
        for r in (1.5, 2.0, math.e, 4.0):
            t_star = (r - 1.0) / math.log(r)
            self.assertAlmostEqual(
                log_gap(2.0, r, t_star),
                span_information(r, 2.0),
                places=12,
            )
            grid = np.linspace(1.0, r, 400)
            numerical = max(abs(log_gap(2.0, r, float(t))) for t in grid)
            self.assertAlmostEqual(numerical, span_information(r, 2.0), delta=1e-4)

    def test_endpoints_vanish(self) -> None:
        self.assertAlmostEqual(log_gap(2.0, 3.0, 1.0), 0.0, places=12)
        self.assertAlmostEqual(log_gap(2.0, 3.0, 3.0), 0.0, places=12)

    def test_small_span_leading_term(self) -> None:
        r = 1.05
        exact = span_information(r, 2.0)
        leading = small_span_leading(r, 2.0)
        self.assertLess(abs(exact - leading) / leading, 0.08)

    def test_span_table_values(self) -> None:
        table = {
            1.5: 0.0205,
            2.0: 0.0597,
            math.e: 0.123,
            3.0: 0.148,
            4.0: 0.234,
            8.0: 0.511,
        }
        for r, expected in table.items():
            self.assertAlmostEqual(span_information(r, 1.0), expected, places=3)

    def test_two_x_critical_n(self) -> None:
        n_star = critical_sample_size(2.0, 2.0)
        self.assertGreater(n_star, 400)
        self.assertLess(n_star, 800)
        self.assertAlmostEqual(
            span_information(2.0, 2.0) * math.sqrt(n_star),
            Z_DETECT,
            places=10,
        )


class TestingBoundTests(unittest.TestCase):
    def test_affinity_between_zero_and_one(self) -> None:
        aff = hellinger_affinity(2.0, 2.0, 200.0)
        self.assertGreater(aff, 0.0)
        self.assertLessEqual(aff, 1.0)

    def test_larger_span_is_easier(self) -> None:
        risk_short = lecam_risk_lower_bound(2.0, 2.0, 200.0)
        risk_long = lecam_risk_lower_bound(4.0, 2.0, 200.0)
        self.assertGreater(risk_short, risk_long)

    def test_two_x_few_points_unmeasurable(self) -> None:
        status = measurability(n_win=80, span_ratio=2.0, shells=8)
        self.assertEqual(status["verdict"], "unmeasurable")

    def test_wide_span_tree_sized_sample_measurable(self) -> None:
        status = measurability(n_win=400, span_ratio=6.0, shells=8)
        self.assertEqual(status["verdict"], "measurable")

    def test_few_shells_unmeasurable_regardless_of_n(self) -> None:
        status = measurability(n_win=10_000, span_ratio=8.0, shells=3)
        self.assertEqual(status["verdict"], "unmeasurable")
        self.assertIn("shells", status["reason"])


class GateSweepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tree = balanced_binary_tree_distances(depth=8)
        cls.grid = integer_grid_distances(side=17)

    def _windowed_call(self, distances: np.ndarray, span: float) -> dict:
        profile = occupancy_profile(distances)
        assert profile is not None
        restricted = restrict_span(*profile, span)
        assert restricted is not None
        return classify_occupancy(*restricted)

    def test_full_span_tree_is_exponential(self) -> None:
        result = analyze_distances(self.tree[0], k_min=6)
        self.assertEqual(result["gate"]["winner"], "exponential")
        self.assertTrue(result["measurable"])
        self.assertEqual(result["growth_class"], "exponential")

    def test_full_span_grid_is_polynomial(self) -> None:
        result = analyze_distances(self.grid[0], k_min=6)
        self.assertEqual(result["gate"]["winner"], "polynomial")
        self.assertTrue(result["measurable"])
        self.assertEqual(result["growth_class"], "polynomial")

    def test_gate_dies_at_two_x_span(self) -> None:
        """The empirical finding: at 2× span the R² gate is not a classifier.

        On matched 2× windows the two hosts no longer separate. That is the
        theorem's shadow, not a software bug.
        """
        tree_short = self._windowed_call(self.tree[0], span=2.0)
        grid_short = self._windowed_call(self.grid[0], span=2.0)
        self.assertLessEqual(tree_short["span_ratio"], 2.01)
        self.assertLessEqual(grid_short["span_ratio"], 2.01)
        tree_full = classify_occupancy(*occupancy_profile(self.tree[0]))
        grid_full = classify_occupancy(*occupancy_profile(self.grid[0]))
        # Noiseless synthetics may still pick the right winner; the
        # signal that does the picking collapses. That is the death.
        self.assertLess(
            tree_short["adjusted_r_squared_margin"],
            0.5 * tree_full["adjusted_r_squared_margin"],
        )
        self.assertLess(
            grid_short["adjusted_r_squared_margin"],
            0.5 * grid_full["adjusted_r_squared_margin"],
        )
        refused_tree = analyze_distances(self.tree[0], span=2.0)
        refused_grid = analyze_distances(self.grid[0], span=2.0)
        self.assertNotEqual(refused_tree["verdict"], "measurable")
        self.assertNotEqual(refused_grid["verdict"], "measurable")
        self.assertIsNone(refused_tree["growth_class"])
        self.assertIsNone(refused_grid["growth_class"])

    def test_wide_window_recovers_both_hosts(self) -> None:
        tree = self._windowed_call(self.tree[0], span=6.0)
        grid = self._windowed_call(self.grid[0], span=6.0)
        self.assertEqual(tree["winner"], "exponential")
        self.assertEqual(grid["winner"], "polynomial")


if __name__ == "__main__":
    unittest.main()
