#!/usr/bin/env python3
"""E5 small-scale: description-length pressure in a trained hierarchy.

Layer IIb rehearsal, labeled as such. This asks whether utilization is
*speakable* with present instruments (growth-class gate + uncertified
occupancy slope). It does not ask whether genomes saturate.

Firewall: β is read from the generator (tag `generator`), never from the
embedding. Packing and defect are read from the learned metric
(tag `representation_metric`). η is refused unless growth class is
exponential; even then the magnitude is uncertified.

Arms (pre-registered together)
------------------------------
1. tree + Poincaré + tree-distance matching (λ = 1)
2. same, λ ∈ {0.1, 1, 10}
3. tree + Poincaré + shuffled targets (no description-length pressure)
4. 2-D grid + Poincaré + distance matching (no hierarchy)
5. tree + Euclidean R² + tree-distance matching (polynomial-exclusion)

Usage
-----
    python3 experiments/e5_trained_hierarchy.py
    python3 experiments/e5_trained_hierarchy.py --quick
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.addressability_meter import (
    measure_ball_growth,
    measure_quartet_defect,
    state_summary,
    validate_distance_matrix,
)
from tools.growth_class_gate import analyze_distances

SCHEMA_VERSION = "e5.small.1"
KAPPA = 1.0
BALL_EPS = 1e-4


def complete_tree(branching: int, depth: int) -> tuple[np.ndarray, np.ndarray]:
    """Hop-distance matrix and depths for a complete rooted b-ary tree."""
    parents = [-1]
    frontier = [0]
    depths = [0]
    for d in range(1, depth + 1):
        nxt = []
        for node in frontier:
            for _ in range(branching):
                parents.append(node)
                depths.append(d)
                nxt.append(len(parents) - 1)
        frontier = nxt
    n = len(parents)
    ancestors: list[dict[int, int]] = []
    for node in range(n):
        chain: dict[int, int] = {}
        current = node
        dist = 0
        while current >= 0:
            chain[current] = dist
            current = parents[current]
            dist += 1
        ancestors.append(chain)
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            common = set(ancestors[i]).intersection(ancestors[j])
            hop = min(ancestors[i][a] + ancestors[j][a] for a in common)
            dist[i, j] = dist[j, i] = hop
    return dist, np.asarray(depths, dtype=np.float64)


def grid_distances(side: int) -> tuple[np.ndarray, np.ndarray]:
    """ℓ¹ grid, origin at a corner. Depth = graph distance from origin."""
    coords = np.array(
        [(i, j) for i in range(side) for j in range(side)], dtype=np.float64
    )
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sum(np.abs(diff), axis=-1)
    depths = coords.sum(axis=1)
    return dist, depths


def to_ball(z: np.ndarray, eps: float = BALL_EPS) -> np.ndarray:
    nrm = np.linalg.norm(z, axis=-1, keepdims=True)
    scale = np.tanh(nrm) * (1.0 - eps)
    return np.where(nrm > 1e-12, z * (scale / np.clip(nrm, 1e-12, None)), z)


def poincare_pairwise(points: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    u = points[:, None, :]
    v = points[None, :, :]
    diff2 = np.sum((u - v) ** 2, axis=-1)
    u2 = np.sum(points**2, axis=-1)
    v2 = u2
    denom = np.clip((1.0 - u2)[:, None] * (1.0 - v2)[None, :], eps, None)
    arg = np.clip(1.0 + 2.0 * diff2 / denom, 1.0 + eps, None)
    return np.arccosh(arg)


def euclidean_pairwise(points: np.ndarray) -> np.ndarray:
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(diff**2, axis=-1) + 1e-18)


def poincare_radius(points: np.ndarray) -> np.ndarray:
    nrm = np.clip(np.linalg.norm(points, axis=-1), 0.0, 1.0 - 1e-6)
    return 2.0 * np.arctanh(nrm)


def polar_init(depths: np.ndarray, host: str, seed: int) -> np.ndarray:
    """Equal-angle shells. A starting chart, not a fitted host."""
    rng = np.random.default_rng(seed)
    z = rng.normal(0.0, 0.02, size=(len(depths), 2))
    by: dict[int, list[int]] = {}
    for i, d in enumerate(depths):
        by.setdefault(int(d), []).append(i)
    for d, nodes in by.items():
        if d == 0:
            z[nodes[0]] = 0.0
            continue
        rho = 0.4 * d
        r = math.tanh(rho / 2.0) if host == "poincare" else 0.35 * d
        jitter = rng.normal(0.0, 0.03, size=len(nodes))
        for k, i in enumerate(nodes):
            theta = 2.0 * math.pi * k / len(nodes) + jitter[k]
            z[i] = [r * math.cos(theta), r * math.sin(theta)]
    return z


def embed(
    target: np.ndarray,
    depths: np.ndarray,
    *,
    host: str,
    lam: float,
    seed: int,
    pin_root: bool = True,
) -> tuple[np.ndarray, float]:
    n = target.shape[0]
    z0 = polar_init(depths, host, seed)
    if pin_root:
        z0[0] = 0.0

    def pack(z: np.ndarray) -> np.ndarray:
        return to_ball(z) if host == "poincare" else z

    def pairwise(p: np.ndarray) -> np.ndarray:
        return poincare_pairwise(p) if host == "poincare" else euclidean_pairwise(p)

    def loss_of(flat: np.ndarray) -> float:
        z = flat.reshape(n, 2).copy()
        if pin_root:
            z[0] = 0.0
        p = pack(z)
        dhat = pairwise(p)
        err = dhat - target
        np.fill_diagonal(err, 0.0)
        return float(lam * np.mean(err**2))

    result = minimize(
        loss_of,
        z0.ravel(),
        method="L-BFGS-B",
        options={"maxiter": 400, "ftol": 1e-10},
    )
    z = result.x.reshape(n, 2)
    if pin_root:
        z[0] = 0.0
    return pack(z), float(result.fun), float(result.fun / lam)


def radial_rate(points: np.ndarray, depths: np.ndarray, host: str) -> float | None:
    mask = depths > 0
    if not np.any(mask):
        return None
    if host == "poincare":
        radii = poincare_radius(points)
    else:
        radii = np.linalg.norm(points, axis=-1)
    ratios = radii[mask] / depths[mask]
    ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
    if ratios.size == 0:
        return None
    return float(np.median(ratios))


def evaluate_embedding(
    points: np.ndarray,
    depths: np.ndarray,
    *,
    host: str,
    beta_nats: float | None,
    seed: int,
) -> dict[str, object]:
    if host == "poincare":
        metric = poincare_pairwise(points)
    else:
        metric = euclidean_pairwise(points)
    np.fill_diagonal(metric, 0.0)
    metric = 0.5 * (metric + metric.T)
    validate_distance_matrix(metric, seed=seed)
    growth = measure_ball_growth(
        metric, centers=min(64, metric.shape[0]), root_index=0, seed=seed
    )
    slope = float(growth["reference_center"]["finite_ball_occupancy_slope"])
    growth_class = analyze_distances(metric[0])
    quartet = measure_quartet_defect(metric, max_quartets=4_000, seed=seed)
    c = radial_rate(points, depths, host)
    winner = growth_class.get("growth_class")
    gate_winner = None
    gate = growth_class.get("gate")
    if isinstance(gate, dict):
        gate_winner = gate.get("winner")
    exponential = bool(
        growth_class.get("measurable") and winner == "exponential"
    )
    host_entropy = slope if exponential else None
    state = state_summary(
        dimension=2.0,
        h_bits=None,
        beta_nats=beta_nats if exponential else None,
        radial_rate=c if exponential else None,
        host_entropy=host_entropy,
        host_entropy_source=(
            "uncertified occupancy slope after exponential growth-class "
            "gate; not independently established packing entropy"
            if exponential
            else None
        ),
        assume_isotropic_hyperbolic=False,
    )
    state["eta_status"] = (
        "uncertified_magnitude"
        if exponential and beta_nats is not None and c is not None
        else "refused"
    )
    return {
        "n": int(metric.shape[0]),
        "host": host,
        "growth_class": growth_class,
        "gate_winner_advisory": gate_winner,
        "occupancy_slope": slope,
        "radial_rate_c": c,
        "quartet_exact_fraction": quartet["exact_fraction"],
        "quartet_delta_q50": quartet["normalized_two_delta_over_s2"]["q50"],
        "state": state,
        "provenance": {
            "beta": "generator" if beta_nats is not None else None,
            "packing": "representation_metric",
        },
    }


def run_arm(
    *,
    name: str,
    target: np.ndarray,
    depths: np.ndarray,
    host: str,
    lam: float,
    seed: int,
    beta_nats: float | None,
    shuffle: bool,
) -> dict[str, object]:
    tgt = target.copy()
    if shuffle:
        rng = np.random.default_rng(seed + 17)
        perm = rng.permutation(tgt.shape[0])
        tgt = tgt[np.ix_(perm, perm)]
    points, loss_weighted, mse = embed(tgt, depths, host=host, lam=lam, seed=seed)
    measured = evaluate_embedding(
        points, depths, host=host, beta_nats=beta_nats, seed=seed
    )
    measured.update(
        {
            "arm": name,
            "lambda": lam,
            "seed": seed,
            "stress": loss_weighted,
            "mse": mse,
            "shuffled_targets": shuffle,
        }
    )
    return measured


def decide(rows: list[dict[str, object]]) -> dict[str, object]:
    def eta(row: dict[str, object]) -> float | None:
        st = row["state"]
        if st.get("eta_status") != "uncertified_magnitude":
            return None
        val = st.get("efficiency_eta")
        return float(val) if val is not None else None

    def winner(row: dict[str, object]) -> str | None:
        official = row["growth_class"].get("growth_class")
        if official in {"exponential", "polynomial", "undecided"}:
            return str(official)
        return None

    def mean_for(arm: str, lam: float | None = None) -> dict[str, object]:
        subset = [r for r in rows if r["arm"] == arm]
        if lam is not None:
            subset = [r for r in subset if r["lambda"] == lam]
        etas = [eta(r) for r in subset]
        etas = [e for e in etas if e is not None]
        wins = [winner(r) for r in subset]
        advisory = [r.get("gate_winner_advisory") for r in subset]
        measurable = [bool(r["growth_class"].get("measurable")) for r in subset]
        return {
            "n": len(subset),
            "measurable_fraction": (
                sum(measurable) / len(measurable) if measurable else 0.0
            ),
            "growth_class_calls": wins,
            "gate_winner_advisory": advisory,
            "eta_uncertified": etas,
            "eta_mean": float(np.mean(etas)) if etas else None,
            "exponential_fraction": (
                sum(w == "exponential" for w in wins) / len(wins) if wins else 0.0
            ),
        }

    arm1 = mean_for("tree_poincare", 1.0)
    arm3 = mean_for("tree_poincare_shuffled", 1.0)
    arm4 = mean_for("grid_poincare", 1.0)
    arm5 = mean_for("tree_euclidean", 1.0)
    sweep = {
        str(lam): mean_for("tree_poincare", lam) for lam in (0.1, 1.0, 10.0)
    }
    etas = [
        sweep[str(lam)]["eta_mean"]
        for lam in (0.1, 1.0, 10.0)
        if sweep[str(lam)]["eta_mean"] is not None
    ]
    sweep_flat = False
    if len(etas) >= 2:
        sweep_flat = (max(etas) - min(etas)) < 0.05

    # Same-data negative control (arm 3) is the load-bearing bake-in test.
    # Arm 4 changes the generator, not only the loss.
    baked_in = (
        arm1["exponential_fraction"] >= 0.5
        and arm3["exponential_fraction"] >= 0.5
    )
    poly_kill = (
        arm5["exponential_fraction"] >= 0.5 and arm5["eta_mean"] is not None
    )
    separated = (
        arm1["exponential_fraction"] >= 0.5
        and arm4["exponential_fraction"] < 0.5
        and arm3["exponential_fraction"] < 0.5
    )
    if baked_in:
        verdict = "KILLED_baked_in_artifact"
        note = (
            "Shuffled-target control (same tree, no correspondence) is still "
            "called exponential. Occupancy is reading the radial tree layout, "
            "not description-length drive. η is not speakable."
        )
    elif poly_kill:
        verdict = "KILLED_polynomial_exclusion"
        note = "Euclidean host appeared to carry exponential retained novelty."
    elif sweep_flat and arm1["exponential_fraction"] >= 0.5:
        verdict = "KILLED_no_rate_response"
        note = (
            "λ sweep is flat. The rate-distortion knob does not move η. "
            "Not a saturation drive."
        )
    elif not separated:
        verdict = "NOT_SPEAKABLE"
        note = (
            "Growth class did not separate hierarchical Poincaré from "
            "controls. Utilization is not a reportable number at this scale."
        )
    else:
        verdict = "SPEAKABLE_rehearsal"
        note = (
            "Arms separate on the growth-class gate. η remains an "
            "uncertified magnitude and is not a saturation result."
        )
    return {
        "verdict": verdict,
        "note": note,
        "lambda_sweep_flat": sweep_flat,
        "arm1_tree_poincare": arm1,
        "arm2_lambda_sweep": sweep,
        "arm3_shuffled": arm3,
        "arm4_grid": arm4,
        "arm5_euclidean_tree": arm5,
        "claim": (
            "IIb instrument rehearsal. Not genomic saturation. "
            "Not a state-equation measurement."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="depth 4, one seed")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/e5_small_scale.json"),
    )
    parser.add_argument(
        "--rescore",
        action="store_true",
        help="recompute the decision from an existing results JSON",
    )
    args = parser.parse_args()

    if args.rescore:
        payload = json.loads(args.output.read_text())
        rows = payload["runs"]
        for row in rows:
            if "mse" not in row and row.get("lambda"):
                row["mse"] = float(row["stress"]) / float(row["lambda"])
        payload["decision"] = decide(rows)
        args.output.write_text(json.dumps(payload, indent=2, default=str))
        print(json.dumps(payload["decision"], indent=2))
        print(f"rescored {args.output}")
        return 0

    depth = 4 if args.quick else 6
    seeds = [0] if args.quick else [0, 1, 2]
    branching = 2
    tree, depths = complete_tree(branching, depth)
    grid, grid_depths = grid_distances(8 if depth >= 6 else 5)
    beta = math.log(branching)

    preregistration = {
        "experiment": "E5_small_scale",
        "layer": "IIb",
        "schema_version": SCHEMA_VERSION,
        "claim": "speakability of utilization, not genomic saturation",
        "tree": {"branching": branching, "depth": depth, "n": int(tree.shape[0])},
        "seeds": seeds,
        "kappa_poincare": KAPPA,
        "lambda_sweep": [0.1, 1.0, 10.0],
        "beta_nats_generator": beta,
        "decision_rule": {
            "predict": (
                "arm1 exponential; arm3 weaker or fail; arm4/euclidean "
                "polynomial or η refused"
            ),
            "kill_baked_in": "arm3 (shuffled targets, same tree) still exponential",
            "kill_polynomial": "euclidean tree arm exponential with reportable η",
            "not_speakable": "growth class does not separate arms",
        },
    }

    rows: list[dict[str, object]] = []
    jobs = []
    for seed in seeds:
        jobs.append(
            ("tree_poincare", tree, depths, "poincare", 0.1, seed, beta, False)
        )
        jobs.append(
            ("tree_poincare", tree, depths, "poincare", 1.0, seed, beta, False)
        )
        jobs.append(
            ("tree_poincare", tree, depths, "poincare", 10.0, seed, beta, False)
        )
        jobs.append(
            (
                "tree_poincare_shuffled",
                tree,
                depths,
                "poincare",
                1.0,
                seed,
                beta,
                True,
            )
        )
        jobs.append(
            (
                "grid_poincare",
                grid,
                grid_depths,
                "poincare",
                1.0,
                seed,
                None,
                False,
            )
        )
        jobs.append(
            (
                "tree_euclidean",
                tree,
                depths,
                "euclidean",
                1.0,
                seed,
                beta,
                False,
            )
        )

    print(f"E5 small-scale: {len(jobs)} embeddings, depth={depth}, seeds={seeds}")
    for i, job in enumerate(jobs, 1):
        name, target, dep, host, lam, seed, beta_nats, shuffle = job
        print(
            f"  [{i}/{len(jobs)}] {name} λ={lam} seed={seed} host={host}",
            flush=True,
        )
        rows.append(
            run_arm(
                name=name,
                target=target,
                depths=dep,
                host=host,
                lam=lam,
                seed=seed,
                beta_nats=beta_nats,
                shuffle=shuffle,
            )
        )

    decision = decide(rows)
    payload = {
        "preregistration": preregistration,
        "runs": rows,
        "decision": decision,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, default=str))
    print(json.dumps(decision, indent=2))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
