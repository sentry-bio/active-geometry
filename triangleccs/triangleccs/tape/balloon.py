"""Yule / JC69 balloon — finite-alphabet tape vs genealogy.

A ground-truth Yule tree with JC69 evolution along branches shows that
endpoint distinguishability (block) can remain high while quartet
resolvability collapses. An event-matched infinite-sites control on the
same branches does not collapse. This animates the ladder; it is not a
biological regime and must not supervise topology.

Sequence simulation walks the tree (shared history). Independent mutation
from a common root with a single p is a star, not a Yule balloon.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from triangleccs.classifier.quartets import (
    classify_quartet_resolvability,
    measure_quartet_defect,
)


@dataclass(frozen=True)
class YuleTree:
    n_tips: int
    root: int
    parent: dict[int, int]
    length: dict[int, float]
    children: dict[int, list[int]]


@dataclass(frozen=True)
class BalloonResult:
    depth: float
    n_tips: int
    seq_len: int
    block_fraction: float
    resolvable_fraction: float
    exact_among_resolvable: float
    model: str


def yule_tree(
    n_tips: int, depth: float, rng: np.random.Generator, birth_rate: float = 1.0
) -> YuleTree:
    """Backward Yule, scaled so every tip is at geodesic depth ``depth`` from the root.

    ``depth`` is expected substitutions on the root→tip path (ultrametric).
    """
    if n_tips < 2:
        raise ValueError("need at least two tips")
    if depth < 0:
        raise ValueError("depth must be non-negative")
    active = list(range(n_tips))
    parent: dict[int, int] = {}
    length: dict[int, float] = {}
    children: dict[int, list[int]] = {i: [] for i in range(n_tips)}
    height_from_present = {i: 0.0 for i in range(n_tips)}
    next_id = n_tips
    time_from_present = 0.0
    while len(active) > 1:
        k = len(active)
        time_from_present += float(rng.exponential(1.0 / (k * birth_rate)))
        a, b = (int(x) for x in rng.choice(active, size=2, replace=False))
        node = next_id
        next_id += 1
        children[node] = [a, b]
        parent[a] = node
        parent[b] = node
        length[a] = time_from_present - height_from_present[a]
        length[b] = time_from_present - height_from_present[b]
        height_from_present[node] = time_from_present
        active.remove(a)
        active.remove(b)
        active.append(node)
    root = int(active[0])
    parent[root] = -1
    length[root] = 0.0
    root_h = height_from_present[root]
    scale = (depth / root_h) if root_h > 0 else 1.0
    length = {k: v * scale for k, v in length.items()}
    return YuleTree(
        n_tips=n_tips, root=root, parent=parent, length=length, children=children
    )


def _path_to_root(tree: YuleTree, tip: int) -> dict[int, float]:
    dist = 0.0
    node = tip
    out: dict[int, float] = {}
    while node != -1:
        out[node] = dist
        dist += tree.length[node]
        node = tree.parent[node]
    return out


def pairwise_tree_distance(tree: YuleTree) -> np.ndarray:
    n = tree.n_tips
    d = np.zeros((n, n), dtype=np.float64)
    paths = [_path_to_root(tree, i) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            node = j
            acc = 0.0
            while node not in paths[i]:
                acc += tree.length[node]
                node = tree.parent[node]
            dist = acc + paths[i][node]
            d[i, j] = d[j, i] = dist
    return d


def simulate_jc69_on_tree(
    tree: YuleTree, seq_len: int, rng: np.random.Generator
) -> np.ndarray:
    """JC69 along branches. Returns tip sequences shape (n_tips, L), bases 0..3."""
    seqs: dict[int, np.ndarray] = {
        tree.root: rng.integers(0, 4, size=seq_len, dtype=np.int8)
    }

    def walk(node: int) -> None:
        for child in tree.children.get(node, []):
            t = tree.length[child]
            p = 0.75 * (1.0 - math.exp(-4.0 * t / 3.0))
            parent_seq = seqs[node]
            mut = rng.random(seq_len) < p
            child_seq = parent_seq.copy()
            n_mut = int(mut.sum())
            if n_mut:
                child_seq[mut] = (child_seq[mut] + rng.integers(1, 4, size=n_mut)) % 4
            seqs[child] = child_seq
            walk(child)

    walk(tree.root)
    return np.stack([seqs[i] for i in range(tree.n_tips)])


def infinite_sites_distance(
    tree: YuleTree, rng: np.random.Generator, seq_len: int = 1
) -> np.ndarray:
    """Event-matched infinite-sites: Poisson unique mutations on each branch.

    Pairwise distance is the mutation count on the path (additive; no ceiling).
    ``seq_len`` scales the Poisson mean so the control is matched to a JC
    alignment of that length (expected events ≈ L · t per branch).
    """
    muts = {tree.root: 0}

    def walk(node: int) -> None:
        for child in tree.children.get(node, []):
            muts[child] = int(rng.poisson(tree.length[child] * seq_len))
            walk(child)

    walk(tree.root)
    n = tree.n_tips
    d = np.zeros((n, n), dtype=np.float64)
    paths = [_path_to_root(tree, i) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            node = j
            acc = 0
            while node not in paths[i]:
                acc += muts[node]
                node = tree.parent[node]
            # mutations from i up to LCA
            acc_i = 0
            climb = i
            while climb != node:
                acc_i += muts[climb]
                climb = tree.parent[climb]
            d[i, j] = d[j, i] = float(acc + acc_i)
    return d


def jc69_p_distance(expected_subs: np.ndarray) -> np.ndarray:
    """JC69 expected p-distance: p = 0.75 (1 - exp(-4D/3))."""
    d = np.asarray(expected_subs, dtype=np.float64)
    p = 0.75 * (1.0 - np.exp(-4.0 * d / 3.0))
    if p.ndim == 2:
        np.fill_diagonal(p, 0.0)
    return p


def jc69_corrected(p: np.ndarray) -> np.ndarray:
    """JC69 distance correction; saturates as p → 0.75."""
    p = np.asarray(p, dtype=np.float64)
    x = np.clip(1.0 - 4.0 * p / 3.0, 1e-12, None)
    d = -0.75 * np.log(x)
    if d.ndim == 2:
        np.fill_diagonal(d, 0.0)
    return d


def hamming_p(tips: np.ndarray) -> np.ndarray:
    tips = np.asarray(tips)
    n, _L = tips.shape
    neq = tips[:, None, :] != tips[None, :, :]
    p = neq.mean(axis=-1).astype(np.float64)
    np.fill_diagonal(p, 0.0)
    return p


def block_fraction_sequences(tips: np.ndarray) -> float:
    """Fraction of unique tip sequences (ε-separation at sequence equality)."""
    n = tips.shape[0]
    uniq = np.unique(tips, axis=0).shape[0]
    return uniq / n


def run_balloon_cell(
    *,
    n_tips: int = 32,
    depth: float = 2.0,
    seq_len: int = 500,
    seed: int = 0,
    max_quartets: int = 800,
    jc_ceiling: float = 0.74,
) -> dict[str, object]:
    """One depth cell: JC69 balloon vs infinite-sites control on the same Yule tree."""
    rng = np.random.default_rng(seed)
    tree = yule_tree(n_tips, depth, rng)
    truth = pairwise_tree_distance(tree)
    tips = simulate_jc69_on_tree(tree, seq_len, rng)
    p = hamming_p(tips)
    d_jc = jc69_corrected(p)
    block = block_fraction_sequences(tips)
    jc_stats = classify_quartet_resolvability(
        p,
        jc_ceiling=jc_ceiling,
        max_quartets=max_quartets,
        seed=seed,
        truth_matrix=truth,
        topology_matrix=d_jc,
    )
    inf = infinite_sites_distance(tree, rng, seq_len=seq_len)
    inf_stats = classify_quartet_resolvability(
        inf,
        jc_ceiling=None,
        max_quartets=max_quartets,
        seed=seed,
        truth_matrix=truth,
        topology_matrix=inf,
    )
    jc = BalloonResult(
        depth=depth,
        n_tips=n_tips,
        seq_len=seq_len,
        block_fraction=block,
        resolvable_fraction=float(jc_stats["resolvable_fraction"]),
        exact_among_resolvable=float(jc_stats["exact_among_resolvable"]),
        model="jc69",
    )
    inf_res = BalloonResult(
        depth=depth,
        n_tips=n_tips,
        seq_len=seq_len,
        block_fraction=1.0,
        resolvable_fraction=float(inf_stats["resolvable_fraction"]),
        exact_among_resolvable=float(inf_stats["exact_among_resolvable"]),
        model="infinite_sites",
    )
    return {"jc69": jc, "infinite_sites": inf_res, "tree": tree, "truth": truth}


def resolution_boundary_depth(
    *,
    seq_len: int,
    n_tips: int = 16,
    threshold: float = 0.5,
    depths: tuple[float, ...] = (0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0),
    seed: int = 0,
    max_quartets: int = 400,
) -> float:
    """Smallest tested depth at which JC resolvable fraction falls below ``threshold``.

    INSTRUMENT. The reported Yule/JC69 sweep found D* ∝ log L; this helper is
    the finite-sample analogue on a coarse grid, not that published fit.
    """
    for depth in depths:
        cell = run_balloon_cell(
            n_tips=n_tips,
            depth=float(depth),
            seq_len=seq_len,
            seed=seed,
            max_quartets=max_quartets,
        )
        if cell["jc69"].resolvable_fraction < threshold:
            return float(depth)
    return float(depths[-1])


def fit_dstar_logL(
    lengths: tuple[int, ...],
    dstar: tuple[float, ...],
) -> dict[str, float | str]:
    """Ordinary least squares D* = a + b ln L. Tests the fitter, not biology."""
    L = np.asarray(lengths, dtype=np.float64)
    y = np.asarray(dstar, dtype=np.float64)
    if L.size != y.size or L.size < 3:
        raise ValueError("need at least three matched (L, D*) points")
    x = np.log(L)
    design = np.column_stack([np.ones(L.size), x])
    a, b = np.linalg.lstsq(design, y, rcond=None)[0]
    pred = a + b * x
    residual = float(np.sum((y - pred) ** 2))
    total = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - residual / total if total > 0 else 1.0
    return {
        "intercept": float(a),
        "slope_lnL": float(b),
        "r_squared": r2,
        "status": "INSTRUMENT",
        "note": "D* = a + b ln L on supplied points; not a clade measurement",
    }


def phase0_certify(seed: int = 0) -> dict[str, object]:
    """Sanity controls: exact Yule tree δ ≈ 0; noise not; JC formula check."""
    rng = np.random.default_rng(seed)
    tree = yule_tree(16, 1.0, rng)
    ultra = pairwise_tree_distance(tree)
    exact = measure_quartet_defect(ultra, max_quartets=500, seed=seed)
    noise = rng.random((16, 16))
    noise = (noise + noise.T) / 2.0
    np.fill_diagonal(noise, 0.0)
    bad = measure_quartet_defect(noise, max_quartets=500, seed=seed)

    depths = np.array([0.05, 0.5, 1.0, 3.0, 5.0])
    p_direct = jc69_p_distance(depths)
    max_dev = float(
        np.max(np.abs(p_direct - 0.75 * (1.0 - np.exp(-4.0 * depths / 3.0))))
    )

    return {
        "exact_tree_exact_fraction": exact["exact_fraction"],
        "noise_exact_fraction": bad["exact_fraction"],
        "jc_formula_max_dev": max_dev,
        "passed": bool(
            float(exact["exact_fraction"]) > 0.95
            and float(bad["exact_fraction"]) < 0.5
            and max_dev < 1e-9
        ),
    }
