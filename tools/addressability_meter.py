#!/usr/bin/env python3
"""Small, encoder-free meter for Active Geometry.

Input
-----
A finite metric distance matrix in ``.npy`` or numeric CSV format.

Output
------
1. A finite-sample ball-occupancy slope, not silently promoted to packing
   entropy.
2. A four-point tree defect, independent of the growth estimate.
3. Optional addressability efficiency when an independently measured process
   rate and radial calibration are supplied.

The program deliberately does not infer a missing process rate or radial
calibration from curvature. It therefore cannot manufacture state-equation
agreement by back-solving one axis from the other.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np

SCHEMA_VERSION = "1.0"


def load_distance_matrix(path: Path) -> np.ndarray:
    """Load a square numeric matrix from NPY or CSV."""
    suffix = path.suffix.lower()
    if suffix == ".npy":
        matrix = np.load(path, allow_pickle=False)
    elif suffix in {".csv", ".txt", ".tsv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        matrix = np.loadtxt(path, delimiter=delimiter)
    else:
        raise ValueError("distance matrix must be .npy, .csv, .tsv, or .txt")
    return np.asarray(matrix, dtype=np.float64)


def validate_distance_matrix(
    matrix: np.ndarray,
    *,
    tolerance: float = 1e-9,
    triangle_samples: int = 10_000,
    seed: int = 0,
) -> dict[str, object]:
    """Validate metric axioms, recording whether triangle checks were sampled."""
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"distance matrix must be square; got {matrix.shape}")
    if matrix.shape[0] < 4:
        raise ValueError("at least four points are required")
    if not np.isfinite(matrix).all():
        raise ValueError("distance matrix contains NaN or infinity")
    scale = float(matrix.max())
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("distance matrix must contain positive distances")
    absolute_tolerance = tolerance * scale
    if float(matrix.min()) < -absolute_tolerance:
        raise ValueError("distance matrix contains negative distances")
    if not np.allclose(
        np.diag(matrix), 0.0, atol=absolute_tolerance, rtol=0.0
    ):
        raise ValueError("distance matrix diagonal must be zero")
    if not np.allclose(
        matrix, matrix.T, atol=absolute_tolerance, rtol=0.0
    ):
        raise ValueError("distance matrix must be symmetric")
    off_diagonal = ~np.eye(matrix.shape[0], dtype=bool)
    if np.any(matrix[off_diagonal] <= 0):
        raise ValueError("distinct points must have strictly positive distance")

    n = matrix.shape[0]
    total = n**3
    if total <= triangle_samples:
        triples: Iterable[tuple[int, int, int]] = itertools.product(
            range(n), repeat=3
        )
        policy = "exhaustive"
        checked = total
    else:
        if triangle_samples <= 0:
            raise ValueError("triangle_samples must be positive")
        rng = np.random.default_rng(seed)
        triples = (
            tuple(int(x) for x in rng.integers(0, n, size=3))
            for _ in range(triangle_samples)
        )
        policy = "sampled"
        checked = triangle_samples
    for i, j, k in triples:
        if (
            matrix[i, k]
            > matrix[i, j] + matrix[j, k] + absolute_tolerance
        ):
            raise ValueError(
                f"triangle inequality violated at ({i}, {j}, {k})"
            )
    return {
        "metric_axioms_checked": True,
        "triangle_policy": policy,
        "triangle_checks": checked,
        "relative_tolerance": tolerance,
    }


def _quantiles(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {"q05": math.nan, "q50": math.nan, "q95": math.nan}
    q05, q50, q95 = np.quantile(values, [0.05, 0.5, 0.95])
    return {"q05": float(q05), "q50": float(q50), "q95": float(q95)}


def _quartets(
    n: int, max_quartets: int, seed: int
) -> Iterator[tuple[int, int, int, int]]:
    if max_quartets <= 0:
        raise ValueError("max_quartets must be positive")
    total = math.comb(n, 4)
    if total <= max_quartets:
        yield from itertools.combinations(range(n), 4)
        return

    rng = np.random.default_rng(seed)
    selected: set[tuple[int, int, int, int]] = set()
    max_attempts = max_quartets * 20
    for _ in range(max_attempts):
        quartet = tuple(sorted(int(x) for x in rng.choice(n, 4, replace=False)))
        selected.add(quartet)
        if len(selected) == max_quartets:
            break
    if len(selected) != max_quartets:
        raise RuntimeError("could not draw the requested number of unique quartets")
    yield from sorted(selected)


def measure_quartet_defect(
    matrix: np.ndarray,
    *,
    max_quartets: int = 100_000,
    seed: int = 0,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Measure the standard Buneman/Gromov four-point defect.

    For each quartet, sort the three four-point sums as s0 <= s1 <= s2.
    Exact additive tree metrics have s2 - s1 = 0.

    The quartet hyperbolicity defect is δ_q = (s2 - s1) / 2. The normalized
    defect 2δ_q/s2 is scale-invariant. Neither estimates curvature magnitude.
    """
    deltas: list[float] = []
    normalized: list[float] = []
    exact = 0

    for a, b, c, d in _quartets(matrix.shape[0], max_quartets, seed):
        sums = np.sort(
            np.array(
                [
                    matrix[a, b] + matrix[c, d],
                    matrix[a, c] + matrix[b, d],
                    matrix[a, d] + matrix[b, c],
                ],
                dtype=np.float64,
            )
        )
        s0, s1, s2 = (float(x) for x in sums)
        gap = max(0.0, s2 - s1)
        delta = gap / 2.0
        deltas.append(delta)
        normalized.append(gap / s2)
        if gap <= tolerance * s2:
            exact += 1

    delta_array = np.asarray(deltas)
    normalized_array = np.asarray(normalized)
    count = int(delta_array.size)
    return {
        "quartets": count,
        "exact_fraction": exact / count,
        "delta_q": _quantiles(delta_array),
        "normalized_two_delta_over_s2": _quantiles(normalized_array),
        "definition": {
            "delta_q": "(largest sum - second-largest sum) / 2",
            "normalized": "2 * delta_q / largest sum",
        },
    }


def fit_ball_growth(
    distances: Sequence[float],
    *,
    min_count: int = 3,
    max_fraction: float = 0.5,
    points: int = 32,
) -> tuple[float, float, int] | None:
    """Fit log ball occupancy against radius for one center.

    The fit excludes the saturated half of a finite sample. It is a finite
    scale diagnostic, not an asymptotic entropy theorem.
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
    scaled_radii = unique_radii / radius_scale
    design = np.column_stack([np.ones(unique_radii.size), scaled_radii])
    response = np.log(occupancies.astype(np.float64))
    intercept, scaled_slope = np.linalg.lstsq(design, response, rcond=None)[0]
    prediction = intercept + scaled_slope * scaled_radii
    slope = scaled_slope / radius_scale
    residual = float(np.sum((response - prediction) ** 2))
    total = float(np.sum((response - response.mean()) ** 2))
    r_squared = 1.0 - residual / total if total > 0 else 1.0
    return float(max(0.0, slope)), r_squared, int(unique_radii.size)


def measure_ball_growth(
    matrix: np.ndarray,
    *,
    centers: int = 64,
    root_index: int | None = None,
    seed: int = 0,
    max_fraction: float = 0.5,
) -> dict[str, object]:
    """Estimate finite-sample packing growth from a reference and other centers.

    Volume entropy is basepoint-independent only in the asymptotic setting.
    Finite trees have severe boundary bias, so the state proxy uses an explicit
    root when supplied and otherwise the metric medoid. Other centers quantify
    sensitivity rather than being silently averaged into the state estimate.
    """
    n = matrix.shape[0]
    if root_index is not None and not 0 <= root_index < n:
        raise ValueError(f"root_index must be in [0, {n})")
    reference_index = (
        int(root_index)
        if root_index is not None
        else int(np.argmin(matrix.mean(axis=1)))
    )
    reference_policy = "supplied root" if root_index is not None else "metric medoid"

    rng = np.random.default_rng(seed)
    if centers <= 0:
        raise ValueError("centers must be positive")
    if centers >= n:
        center_indices = np.arange(n)
    else:
        candidates = np.delete(np.arange(n), reference_index)
        additional = rng.choice(
            candidates, max(0, centers - 1), replace=False
        )
        center_indices = np.sort(
            np.concatenate([[reference_index], additional])
        )

    slopes: list[float] = []
    fits: list[float] = []
    points: list[int] = []
    used: list[int] = []
    for center in center_indices:
        fitted = fit_ball_growth(
            matrix[center],
            max_fraction=max_fraction,
        )
        if fitted is None:
            continue
        slope, r_squared, fitted_points = fitted
        slopes.append(slope)
        fits.append(r_squared)
        points.append(fitted_points)
        used.append(int(center))

    if not slopes:
        raise ValueError("insufficient distinct radii to estimate ball growth")
    if reference_index not in used:
        raise ValueError("reference center has insufficient radii for growth fit")

    slope_array = np.asarray(slopes)
    fit_array = np.asarray(fits)
    reference_position = used.index(reference_index)
    return {
        "centers_requested": int(center_indices.size),
        "centers_fitted": len(slopes),
        "center_indices": used,
        "reference_center": {
            "index": reference_index,
            "policy": reference_policy,
            "finite_ball_occupancy_slope": slopes[reference_position],
            "fit_r_squared": fits[reference_position],
            "fit_points": points[reference_position],
        },
        "finite_ball_occupancy_slope": _quantiles(slope_array),
        "fit_r_squared": _quantiles(fit_array),
        "median_fit_points": float(np.median(points)),
        "fit_window": {
            "minimum_ball_count": 3,
            "maximum_occupancy_fraction": max_fraction,
        },
        "warning": (
            "finite-sample slope of observed log ball occupancy; this is not "
            "host packing entropy without a separation-scale and depth/size "
            "convergence study"
        ),
    }


def state_summary(
    *,
    dimension: float,
    h_bits: float | None,
    beta_nats: float | None,
    radial_rate: float | None,
    host_entropy: float | None,
    host_entropy_source: str | None,
    assume_isotropic_hyperbolic: bool,
) -> dict[str, object]:
    """Calculate only quantities licensed by independently supplied inputs."""
    values = {
        "dimension": dimension,
        "h_bits": h_bits,
        "beta_nats": beta_nats,
        "radial_rate": radial_rate,
        "host_entropy": host_entropy,
    }
    for name, value in values.items():
        if value is not None and not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if assume_isotropic_hyperbolic and dimension <= 1:
        raise ValueError("isotropic hyperbolic dimension must be greater than one")
    if h_bits is not None and beta_nats is not None:
        raise ValueError("supply h_bits or beta_nats, not both")
    if h_bits is not None and h_bits <= 0:
        raise ValueError("h_bits must be positive")
    if beta_nats is not None and beta_nats <= 0:
        raise ValueError("beta_nats must be positive")
    if radial_rate is not None and radial_rate <= 0:
        raise ValueError("radial_rate must be positive")
    if host_entropy is not None and host_entropy <= 0:
        raise ValueError("host_entropy must be positive")

    beta = beta_nats
    if h_bits is not None:
        beta = h_bits * math.log(2.0)

    result: dict[str, object] = {"dimension": dimension}
    if beta is not None:
        result["beta_nats_per_step"] = beta
    if host_entropy is not None:
        result["host_entropy_nats_per_distance"] = host_entropy
        result["host_entropy_source"] = host_entropy_source
    if host_entropy is not None and radial_rate is not None:
        capacity = radial_rate * host_entropy
        result["capacity_nats_per_step"] = capacity
    if beta is not None and radial_rate is not None and host_entropy is not None:
        capacity = radial_rate * host_entropy
        result.update(
            {
                "efficiency_eta": beta / capacity,
                "addressability_slack": capacity - beta,
                "bound_satisfied": beta <= capacity,
                "equality_axes_supplied": True,
            }
        )
    else:
        result["equality_axes_supplied"] = False
        result["missing_for_equality_test"] = [
            name
            for name, value in (
                ("process_rate", beta),
                ("radial_rate_c", radial_rate),
                ("host_entropy", host_entropy),
            )
            if value is None
        ]

    result["isotropic_hyperbolic_assumption"] = assume_isotropic_hyperbolic
    if assume_isotropic_hyperbolic:
        if host_entropy is not None:
            curvature = (host_entropy / (dimension - 1.0)) ** 2
            result["equivalent_curvature_magnitude"] = curvature
            if radial_rate is not None:
                result["equivalent_normalized_curvature"] = (
                    radial_rate**2 * curvature
                )
        if beta is not None and radial_rate is not None:
            result["isotropic_curvature_floor"] = (
                beta / (radial_rate * (dimension - 1.0))
            ) ** 2
    return result


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def analyze(
    path: Path,
    *,
    dimension: float = 2.0,
    h_bits: float | None = None,
    beta_nats: float | None = None,
    radial_rate: float | None = None,
    host_entropy: float | None = None,
    promote_occupancy_slope: bool = False,
    assume_isotropic_hyperbolic: bool = False,
    max_quartets: int = 100_000,
    centers: int = 64,
    root_index: int | None = None,
    seed: int = 0,
) -> dict[str, object]:
    matrix = load_distance_matrix(path)
    validation = validate_distance_matrix(matrix, seed=seed)
    growth = measure_ball_growth(
        matrix,
        centers=centers,
        root_index=root_index,
        seed=seed,
    )
    occupancy_slope = float(
        growth["reference_center"][  # type: ignore[index]
            "finite_ball_occupancy_slope"
        ]
    )
    if host_entropy is not None and promote_occupancy_slope:
        raise ValueError(
            "supply host_entropy or promote the occupancy slope, not both"
        )
    entropy_for_state = occupancy_slope if promote_occupancy_slope else host_entropy
    entropy_source = None
    if host_entropy is not None:
        entropy_source = "explicit independent input"
    elif promote_occupancy_slope:
        entropy_source = (
            "finite occupancy slope promoted by explicit user assumption; "
            "not independently established packing entropy"
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "input": {
            "path": str(path),
            "sha256": file_sha256(path),
            "points": int(matrix.shape[0]),
            "distance_units": "as supplied; no hidden normalization",
            "validation": validation,
        },
        "growth": growth,
        "quartet": measure_quartet_defect(
            matrix,
            max_quartets=max_quartets,
            seed=seed,
        ),
        "state": state_summary(
            dimension=dimension,
            h_bits=h_bits,
            beta_nats=beta_nats,
            radial_rate=radial_rate,
            host_entropy=entropy_for_state,
            host_entropy_source=entropy_source,
            assume_isotropic_hyperbolic=assume_isotropic_hyperbolic,
        ),
        "independence": {
            "verified": False,
            "note": (
                "The meter records supplied axes but cannot prove that process "
                "rate, radial calibration, and distances were estimated "
                "independently."
            ),
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Measure finite ball growth and four-point tree defect from one "
            "distance matrix."
        )
    )
    parser.add_argument("matrix", type=Path, help=".npy or numeric CSV matrix")
    parser.add_argument("--dimension", type=float, default=2.0)
    rate = parser.add_mutually_exclusive_group()
    rate.add_argument("--h-bits", type=float, help="independent bits/step estimate")
    rate.add_argument(
        "--beta-nats", type=float, help="independent retained growth in nats/step"
    )
    parser.add_argument(
        "--radial-rate",
        type=float,
        help="independent radial distance per generative step (c)",
    )
    entropy = parser.add_mutually_exclusive_group()
    entropy.add_argument(
        "--host-entropy",
        type=float,
        help="independent host packing/volume entropy in nats/distance",
    )
    entropy.add_argument(
        "--promote-occupancy-slope",
        action="store_true",
        help=(
            "explicitly treat the finite occupancy slope as host entropy; "
            "requires an external convergence justification"
        ),
    )
    parser.add_argument(
        "--assume-isotropic-hyperbolic",
        action="store_true",
        help="permit conversion of host entropy to an equivalent curvature",
    )
    parser.add_argument("--max-quartets", type=int, default=100_000)
    parser.add_argument("--centers", type=int, default=64)
    parser.add_argument(
        "--root-index",
        type=int,
        help="root/reference row; default is the metric medoid",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, help="write JSON to this path")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output is not None and args.output.resolve() == args.matrix.resolve():
        raise ValueError("output path must not overwrite the input matrix")
    report = analyze(
        args.matrix,
        dimension=args.dimension,
        h_bits=args.h_bits,
        beta_nats=args.beta_nats,
        radial_rate=args.radial_rate,
        host_entropy=args.host_entropy,
        promote_occupancy_slope=args.promote_occupancy_slope,
        assume_isotropic_hyperbolic=args.assume_isotropic_hyperbolic,
        max_quartets=args.max_quartets,
        centers=args.centers,
        root_index=args.root_index,
        seed=args.seed,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.output:
        args.output.write_text(rendered + "\n")
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
