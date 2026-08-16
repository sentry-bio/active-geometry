"""Sextant v0 — the encoder that belongs to this datum.

A genomic language model compresses strings. The MDL of diversity is a tree
(or its quartets) plus a packing of addresses. The sextant is the map that
takes pairwise distances onto the frozen polar chart. It is not a neural net,
does not set κ, and must not supervise topology from a sequence metric.

See docs/ENCODER.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from triangleccs.address import Address, make_address
from triangleccs.chart.poincare import distance_matrix, exp_map_zero, logmap0
from triangleccs.chart.polar import hyperbolic_radius
from triangleccs.classifier.quartets import classify_quartet_resolvability, measure_quartet_defect
from triangleccs.datum.form import Form
from triangleccs.datum.gauge import wrap_pi
from triangleccs.packing.bound import chart_block_separation
from triangleccs.tape import balloon


BASES = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "N": 0}


def encode_bases(sequences: list[str] | np.ndarray) -> np.ndarray:
    """Aligned DNA/RNA strings or integer array → (n, L) int8 in {0,1,2,3}."""
    if isinstance(sequences, np.ndarray):
        return np.asarray(sequences, dtype=np.int8)
    rows = []
    width = max(len(s) for s in sequences)
    for s in sequences:
        row = np.zeros(width, dtype=np.int8)
        for i, ch in enumerate(s.upper()):
            row[i] = BASES.get(ch, 0)
        rows.append(row)
    return np.stack(rows)


def _cosh_law_angle(ri: float, rj: float, dij: float, kappa: float) -> float:
    sk = math.sqrt(kappa)
    den = math.sinh(sk * ri) * math.sinh(sk * rj)
    if den < 1e-15 or ri <= 0 or rj <= 0:
        return 0.0
    c = (math.cosh(sk * ri) * math.cosh(sk * rj) - math.cosh(sk * dij)) / den
    return float(math.acos(min(1.0, max(-1.0, c))))


def embed_distances_h2(
    distances: np.ndarray,
    form: Form | None = None,
    *,
    origin_index: int = 0,
    meridian_index: int = 1,
) -> np.ndarray:
    """Place points in the Poincaré disk by hyperbolic law of cosines.

    Origin index sits at the chart origin. Meridian index is placed on +x.
    Later points get an angle from the origin–meridian triangle; the sign is
    chosen to match distance to the first placed off-axis point. Residual of
    the chart vs the input metric is the honesty of this embedding, not a
    proof that biology fills H².
    """
    form = form or Form()
    D = np.asarray(distances, dtype=np.float64)
    n = D.shape[0]
    if n < 2:
        raise ValueError("need at least two points")
    r = D[origin_index].copy()
    r[origin_index] = 0.0
    theta = np.zeros(n, dtype=np.float64)
    if meridian_index == origin_index:
        meridian_index = 0 if origin_index != 0 else 1
    theta[meridian_index] = 0.0
    ref: int | None = None
    for i in range(n):
        if i in (origin_index, meridian_index):
            continue
        alpha = _cosh_law_angle(
            float(r[i]), float(r[meridian_index]), float(D[i, meridian_index]), form.kappa
        )
        if ref is None:
            theta[i] = alpha
            ref = i
            continue
        plus = alpha
        minus = -alpha
        d_plus = _chord_geodesic(r[i], r[ref], plus - theta[ref], form.kappa)
        d_minus = _chord_geodesic(r[i], r[ref], minus - theta[ref], form.kappa)
        theta[i] = plus if abs(d_plus - D[i, ref]) <= abs(d_minus - D[i, ref]) else minus
    rho = np.tanh(math.sqrt(form.kappa) * r / 2.0) / math.sqrt(form.kappa)
    xy = np.column_stack([rho * np.cos(theta), rho * np.sin(theta)])
    xy[origin_index] = 0.0
    return xy


def _chord_geodesic(ri: float, rj: float, dtheta: float, kappa: float) -> float:
    sk = math.sqrt(kappa)
    c = math.cosh(sk * ri) * math.cosh(sk * rj) - math.sinh(sk * ri) * math.sinh(
        sk * rj
    ) * math.cos(dtheta)
    c = max(1.0, c)
    return float(math.acosh(c) / sk)


def apply_two_anchor_gauge(
    xy: np.ndarray,
    meridian_index: int,
    chirality_index: int,
    form: Form | None = None,
) -> np.ndarray:
    """Rotate/reflect so meridian is +x and chirality has θ ≥ 0."""
    form = form or Form()
    T = logmap0(xy, form.kappa)
    th = np.arctan2(T[:, 1], T[:, 0])
    th = wrap_pi(th - th[meridian_index])
    if th[chirality_index] < 0:
        th = -th
    r = np.linalg.norm(T, axis=1)
    v = np.column_stack([r * np.cos(th), r * np.sin(th)])
    return exp_map_zero(v, form.kappa)


def chart_residual(xy: np.ndarray, distances: np.ndarray, form: Form | None = None) -> float:
    """RMSE of Poincaré chart distances vs the input metric (INSTRUMENT)."""
    form = form or Form()
    chart = distance_matrix(xy, form)
    D = np.asarray(distances, dtype=np.float64)
    n = D.shape[0]
    if n < 2:
        return 0.0
    iu = np.triu_indices(n, k=1)
    err = chart[iu] - D[iu]
    return float(np.sqrt(np.mean(err * err)))


def per_point_residual(
    xy: np.ndarray, distances: np.ndarray, form: Form | None = None
) -> np.ndarray:
    form = form or Form()
    chart = distance_matrix(xy, form)
    D = np.asarray(distances, dtype=np.float64)
    n = D.shape[0]
    out = np.zeros(n, dtype=np.float64)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        err = chart[i, mask] - D[i, mask]
        out[i] = float(np.sqrt(np.mean(err * err))) if n > 1 else 0.0
    return out


@dataclass(frozen=True)
class SextantReport:
    addresses: list[Address]
    xy: np.ndarray
    p_distance: np.ndarray
    jc_distance: np.ndarray
    delta: float
    resolvable: float
    block_sep: float
    residual: float
    note: str


def place_sequences(
    sequences: list[str] | np.ndarray,
    form: Form | None = None,
    *,
    meridian_index: int = 0,
    chirality_index: int = 1,
    origin_index: int = 0,
    theta_status: str = "candidate",
    max_quartets: int = 2_000,
    seed: int = 0,
) -> SextantReport:
    """Distance sextant: aligned sequences → polar H² addresses with instruments.

    Sequence distances place points; they do not certify topology. Quartet
    fields are INSTRUMENT readings of the received metric.
    """
    form = form or Form()
    tips = encode_bases(sequences)
    if tips.shape[0] < 2:
        raise ValueError("need at least two sequences")
    if meridian_index == chirality_index:
        raise ValueError("meridian and chirality indices must differ")
    p = balloon.hamming_p(tips)
    jc = balloon.jc69_corrected(p)
    xy = embed_distances_h2(
        jc, form, origin_index=origin_index, meridian_index=meridian_index
    )
    # Chart origin is the tangent mean, never a taxon and never LUCA.
    tangent = logmap0(xy, form.kappa)
    xy = exp_map_zero(tangent - tangent.mean(axis=0), form.kappa)
    xy = apply_two_anchor_gauge(xy, meridian_index, chirality_index, form)
    residual = chart_residual(xy, jc, form)
    point_res = per_point_residual(xy, jc, form)
    block = float(chart_block_separation(xy, form))
    if tips.shape[0] >= 4:
        defect = measure_quartet_defect(jc, max_quartets=max_quartets, seed=seed)
        reso = classify_quartet_resolvability(
            p,
            jc_ceiling=0.74,
            max_quartets=max_quartets,
            seed=seed,
            topology_matrix=jc,
        )
        delta = float(defect["delta_q"]["q50"])
        resolvable = float(reso["resolvable_fraction"])
    else:
        delta = 0.0
        resolvable = 1.0
    r = hyperbolic_radius(xy, form)
    tangent = logmap0(xy, form.kappa)
    th = np.arctan2(tangent[:, 1], tangent[:, 0])
    addresses = [
        make_address(
            form=form,
            theta=float(th[i]),
            r=float(r[i]),
            theta_status=theta_status,  # type: ignore[arg-type]
            delta=delta,
            resolvable=resolvable,
            block_sep=block,
            residual=float(point_res[i]),
        )
        for i in range(tips.shape[0])
    ]
    return SextantReport(
        addresses=addresses,
        xy=xy,
        p_distance=p,
        jc_distance=jc,
        delta=delta,
        resolvable=resolvable,
        block_sep=block,
        residual=residual,
        note=(
            "Sextant v0: JC distances onto polar H². Sequence metric must not "
            "supervise topology. r ADVISORY; θ CANDIDATE until freeze-gate; "
            "κ CONVENTION."
        ),
    )
