"""Four-point classifier — tree-ness, not curvature.

Quartets ask whether relations still form a tree. They do not calibrate κ,
capacity, or saturation. Unresolved (channel-exhausted) ≠ confidently wrong.
"""

from __future__ import annotations

import itertools
import math
from typing import Iterator

import numpy as np


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
    for _ in range(max_quartets * 20):
        q = tuple(sorted(int(x) for x in rng.choice(n, 4, replace=False)))
        selected.add(q)
        if len(selected) == max_quartets:
            break
    if len(selected) != max_quartets:
        raise RuntimeError("could not draw requested unique quartets")
    yield from sorted(selected)


def quartet_sums(
    matrix: np.ndarray, a: int, b: int, c: int, d: int
) -> np.ndarray:
    return np.sort(
        np.array(
            [
                matrix[a, b] + matrix[c, d],
                matrix[a, c] + matrix[b, d],
                matrix[a, d] + matrix[b, c],
            ],
            dtype=np.float64,
        )
    )


def measure_quartet_defect(
    matrix: np.ndarray,
    *,
    max_quartets: int = 50_000,
    seed: int = 0,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Buneman / Gromov four-point defect. Does not estimate curvature."""
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("distance matrix must be square")
    if matrix.shape[0] < 4:
        raise ValueError("need at least four points")

    deltas: list[float] = []
    normalized: list[float] = []
    exact = 0
    for a, b, c, d in _quartets(matrix.shape[0], max_quartets, seed):
        s0, s1, s2 = (float(x) for x in quartet_sums(matrix, a, b, c, d))
        gap = max(0.0, s2 - s1)
        delta = gap / 2.0
        deltas.append(delta)
        normalized.append(0.0 if s2 <= 0 else gap / s2)
        if gap <= tolerance * max(s2, 1e-15):
            exact += 1

    arr = np.asarray(deltas)
    narr = np.asarray(normalized)
    count = int(arr.size)

    def _q(x: np.ndarray) -> dict[str, float]:
        if x.size == 0:
            return {"q05": math.nan, "q50": math.nan, "q95": math.nan}
        q05, q50, q95 = np.quantile(x, [0.05, 0.5, 0.95])
        return {"q05": float(q05), "q50": float(q50), "q95": float(q95)}

    return {
        "quartets": count,
        "exact_fraction": exact / count if count else math.nan,
        "delta_q": _q(arr),
        "normalized_two_delta_over_s2": _q(narr),
        "status": "INSTRUMENT",
        "note": "classifies tree-ness; does not calibrate kappa or capacity",
    }


def quartet_split(
    matrix: np.ndarray, a: int, b: int, c: int, d: int
) -> frozenset[frozenset[int]]:
    """Four-point inferred split: pairing with the smallest pair-sum."""
    scores = {
        frozenset({frozenset({a, b}), frozenset({c, d})}): float(
            matrix[a, b] + matrix[c, d]
        ),
        frozenset({frozenset({a, c}), frozenset({b, d})}): float(
            matrix[a, c] + matrix[b, d]
        ),
        frozenset({frozenset({a, d}), frozenset({b, c})}): float(
            matrix[a, d] + matrix[b, c]
        ),
    }
    return min(scores, key=scores.get)


def classify_quartet_resolvability(
    matrix: np.ndarray,
    *,
    jc_ceiling: float | None = None,
    max_quartets: int = 5_000,
    seed: int = 0,
    tolerance: float = 1e-9,
    truth_matrix: np.ndarray | None = None,
    topology_matrix: np.ndarray | None = None,
) -> dict[str, object]:
    """Split unresolved (channel-exhausted) from wrong vs correct topology.

    If jc_ceiling is set (e.g. 0.75 for JC p-distance), any pair reaching the
    ceiling marks the quartet unresolved. Otherwise all quartets are treated
    as resolvable.

    If ``truth_matrix`` is given, ``exact_among_resolvable`` is topology
    accuracy against that additive tree, not four-point defect of ``matrix``.
    ``topology_matrix`` is the matrix used to infer the observed split
    (defaults to ``matrix``; pass JC-corrected distances while ceilinging
    on p-distances).
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    topo = matrix if topology_matrix is None else np.asarray(topology_matrix)
    truth = None if truth_matrix is None else np.asarray(truth_matrix, dtype=np.float64)
    n = matrix.shape[0]
    resolvable = 0
    unresolved = 0
    exact = 0
    for a, b, c, d in _quartets(n, max_quartets, seed):
        idxs = (a, b, c, d)
        if jc_ceiling is not None:
            pairs = [matrix[i, j] for i, j in itertools.combinations(idxs, 2)]
            if any(float(p) >= jc_ceiling - 1e-12 for p in pairs):
                unresolved += 1
                continue
        resolvable += 1
        if truth is not None:
            if quartet_split(topo, a, b, c, d) == quartet_split(truth, a, b, c, d):
                exact += 1
        else:
            s0, s1, s2 = (float(x) for x in quartet_sums(topo, a, b, c, d))
            gap = max(0.0, s2 - s1)
            if gap <= tolerance * max(s2, 1e-15):
                exact += 1

    total = resolvable + unresolved
    return {
        "quartets": total,
        "resolvable": resolvable,
        "unresolved": unresolved,
        "resolvable_fraction": resolvable / total if total else math.nan,
        "exact_among_resolvable": exact / resolvable if resolvable else math.nan,
        "status": "INSTRUMENT",
        "note": "unresolved != wrong; channel exhaustion destroys resolvability",
    }


def path_tree_distance(n_leaves: int, branch: float = 1.0) -> np.ndarray:
    """Additive distances on a balanced caterpillar-ish path tree of n leaves.

    Used only in tests: exact tree ⇒ δ = 0.
    """
    # Star of n leaves with unit edge to a shared root projected to leaves:
    # d(i,j) = 2 * branch for i != j.
    d = np.full((n_leaves, n_leaves), 2.0 * branch, dtype=np.float64)
    np.fill_diagonal(d, 0.0)
    return d
