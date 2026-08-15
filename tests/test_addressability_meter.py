import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tools.addressability_meter import (
    analyze,
    fit_ball_growth,
    measure_ball_growth,
    measure_quartet_defect,
    state_summary,
    validate_distance_matrix,
)


def balanced_binary_tree_distances(depth: int) -> np.ndarray:
    """All-node graph distances for a complete binary tree."""
    parents = [-1]
    frontier = [0]
    for _ in range(depth):
        next_frontier = []
        for node in frontier:
            for _ in range(2):
                parents.append(node)
                next_frontier.append(len(parents) - 1)
        frontier = next_frontier

    ancestors: list[dict[int, int]] = []
    for node in range(len(parents)):
        chain: dict[int, int] = {}
        current = node
        distance = 0
        while current >= 0:
            chain[current] = distance
            current = parents[current]
            distance += 1
        ancestors.append(chain)

    matrix = np.zeros((len(parents), len(parents)), dtype=float)
    for i in range(len(parents)):
        for j in range(i + 1, len(parents)):
            common = set(ancestors[i]).intersection(ancestors[j])
            distance = min(ancestors[i][a] + ancestors[j][a] for a in common)
            matrix[i, j] = matrix[j, i] = distance
    return matrix


class AddressabilityMeterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tree = balanced_binary_tree_distances(depth=7)

    def test_exact_tree_has_zero_four_point_defect(self) -> None:
        result = measure_quartet_defect(
            self.tree, max_quartets=5_000, seed=3
        )
        self.assertEqual(result["exact_fraction"], 1.0)
        self.assertEqual(result["delta_q"]["q95"], 0.0)
        self.assertEqual(result["normalized_two_delta_over_s2"]["q95"], 0.0)

    def test_non_tree_metric_has_positive_defect(self) -> None:
        points = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
        matrix = np.linalg.norm(points[:, None] - points[None, :], axis=-1)
        validate_distance_matrix(matrix)
        result = measure_quartet_defect(matrix)
        self.assertEqual(result["exact_fraction"], 0.0)
        self.assertGreater(result["normalized_two_delta_over_s2"]["q50"], 0.0)

    def test_binary_tree_root_growth_approaches_log_two(self) -> None:
        fitted = fit_ball_growth(
            self.tree[0], min_count=3, max_fraction=0.5, points=32
        )
        self.assertIsNotNone(fitted)
        slope, r_squared, _ = fitted  # type: ignore[misc]
        self.assertAlmostEqual(slope, math.log(2), delta=0.12)
        self.assertGreater(r_squared, 0.98)

    def test_radial_rescaling_preserves_normalized_outputs(self) -> None:
        growth = measure_ball_growth(self.tree, centers=16, seed=7)
        h1 = growth["finite_ball_occupancy_slope"]["q50"]
        first = state_summary(
            dimension=2,
            h_bits=1.0,
            beta_nats=None,
            radial_rate=1.0,
            host_entropy=h1,
            host_entropy_source="test",
            assume_isotropic_hyperbolic=True,
        )
        for factor in (3.0, 1e-100, 1e100):
            scaled_matrix = self.tree * factor
            validate_distance_matrix(scaled_matrix)
            scaled_growth = measure_ball_growth(
                scaled_matrix, centers=16, seed=7
            )
            h2 = scaled_growth["finite_ball_occupancy_slope"]["q50"]
            self.assertAlmostEqual(h1, factor * h2, places=10)
            self.assertEqual(
                measure_quartet_defect(
                    scaled_matrix, max_quartets=1_000, seed=5
                )["exact_fraction"],
                1.0,
            )

            second = state_summary(
                dimension=2,
                h_bits=1.0,
                beta_nats=None,
                radial_rate=factor,
                host_entropy=h2,
                host_entropy_source="test",
                assume_isotropic_hyperbolic=True,
            )
            self.assertAlmostEqual(
                first["equivalent_normalized_curvature"],
                second["equivalent_normalized_curvature"],
                places=10,
            )
            self.assertAlmostEqual(
                first["efficiency_eta"], second["efficiency_eta"], places=10
            )

    def test_no_rate_means_no_equality_claim(self) -> None:
        result = state_summary(
            dimension=2,
            h_bits=None,
            beta_nats=None,
            radial_rate=None,
            host_entropy=None,
            host_entropy_source=None,
            assume_isotropic_hyperbolic=False,
        )
        self.assertFalse(result["equality_axes_supplied"])
        self.assertNotIn("efficiency_eta", result)
        self.assertNotIn("equivalent_curvature_magnitude", result)

    def test_end_to_end_report_records_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tree.npy"
            np.save(path, self.tree)
            report = analyze(
                path,
                h_bits=1.0,
                radial_rate=1.0,
                max_quartets=1_000,
                centers=8,
                root_index=0,
                promote_occupancy_slope=True,
                assume_isotropic_hyperbolic=True,
                seed=11,
            )
        self.assertEqual(report["schema_version"], "1.0")
        self.assertEqual(report["input"]["points"], self.tree.shape[0])
        self.assertEqual(len(report["input"]["sha256"]), 64)
        self.assertFalse(report["independence"]["verified"])
        self.assertEqual(report["quartet"]["exact_fraction"], 1.0)
        self.assertEqual(report["growth"]["reference_center"]["policy"], "supplied root")
        self.assertTrue(report["state"]["bound_satisfied"])
        self.assertIn("explicit user assumption", report["state"]["host_entropy_source"])
        self.assertIn("growth_class", report)
        self.assertIn(report["growth_class"]["verdict"],
                      ("measurable", "undecided", "unmeasurable"))


if __name__ == "__main__":
    unittest.main()
