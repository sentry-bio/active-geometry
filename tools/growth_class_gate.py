#!/usr/bin/env python3
"""Layer 0 growth-class measurability and classification.

A finite pointed sample can be asked whether its occupancy is exponential
or polynomial only when the window carries enough span information. This
module computes that predicate from the identities in
``theory/MEASURABILITY.md`` and, when the sample is measurable, runs the
adjusted-R² gate.

The module refuses to emit a growth class on unmeasurable data. That
refusal is the theorem's operational content.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

# Defaults match theory/MEASURABILITY.md §3.
DEFAULT_D_STAR = 2.0
DEFAULT_ALPHA = 0.05
DEFAULT_K_MIN = 6
# Composite-test allowance used for the N table in the note:
# demand Δ√N ≥ Z_DETECT before calling a cell "measurable".
Z_DETECT = 3.0


def midpoint_exponent(r: float, d: float) -> float:
    """Lemma 1.2: f(√r) = d log r (1/2 - 1/(√r+1))."""
    if r <= 1.0:
        raise ValueError("radial span ratio r must exceed 1")
    if d <= 0.0:
        raise ValueError("dimension d must be positive")
    return d * math.log(r) * (0.5 - 1.0 / (math.sqrt(r) + 1.0))


def midpoint_ratio(r: float, d: float) -> float:
    """γ(r,d) = exp(f(√r))."""
    return math.exp(midpoint_exponent(r, d))


def span_information(r: float, d: float) -> float:
    """Proposition 1.3: exact maximum log-gap Δ(r,d)."""
    if r <= 1.0:
        raise ValueError("radial span ratio r must exceed 1")
    if d <= 0.0:
        raise ValueError("dimension d must be positive")
    log_r = math.log(r)
    return d * (math.log((r - 1.0) / log_r) - 1.0 + log_r / (r - 1.0))


def small_span_leading(r: float, d: float) -> float:
    """Leading term of Proposition 1.4: d(r-1)²/8."""
    return d * (r - 1.0) ** 2 / 8.0


def critical_sample_size(r: float, d: float, z: float = Z_DETECT) -> float:
    """N such that Δ√N = z. The scaling table in the note."""
    delta = span_information(r, d)
    if delta <= 0.0:
        return math.inf
    return (z / delta) ** 2


def occupancy_profile(
    distances: Sequence[float],
    *,
    min_count: int = 3,
    max_fraction: float = 0.5,
    points: int = 32,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Radii and occupancies used by the exponential/polynomial fits.

    The window excludes the saturated half of a finite sample, matching
    ``fit_ball_growth`` in ``tools/addressability_meter.py``.
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
    return unique_radii.astype(np.float64), occupancies.astype(np.float64)


def restrict_span(
    radii: np.ndarray,
    occupancies: np.ndarray,
    span: float,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Keep the outermost contiguous window whose radial ratio is at most ``span``.

    Starting from the innermost radius of a unit-edge tree yields only two
    shells at span 2 (radii 1 and 2). The outermost window is the one a
    finite sample actually has when its observed range is a factor of
    ``span``.
    """
    if radii.size < 3 or span <= 1.0:
        return None
    outer = float(radii[-1])
    if outer <= 0.0:
        return None
    keep = radii >= outer / span
    if int(np.count_nonzero(keep)) < 3:
        return None
    return radii[keep], occupancies[keep]


def _adjusted_r2(
    response: np.ndarray, prediction: np.ndarray, n_params: int
) -> float:
    n = int(response.size)
    residual = float(np.sum((response - prediction) ** 2))
    total = float(np.sum((response - response.mean()) ** 2))
    if total <= 0.0:
        return 1.0
    r2 = 1.0 - residual / total
    denom = n - n_params - 1
    if denom <= 0:
        return math.nan
    return 1.0 - (1.0 - r2) * (n - 1) / denom


def _fit_affine(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    design = np.column_stack([np.ones(x.size), x])
    intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
    prediction = intercept + slope * x
    return float(intercept), float(slope), _adjusted_r2(y, prediction, 2)


def classify_occupancy(
    radii: np.ndarray, occupancies: np.ndarray
) -> dict[str, object]:
    """Adjusted-R² comparison of log n ~ ρ versus log n ~ log ρ."""
    if radii.size != occupancies.size or radii.size < 3:
        raise ValueError("radii and occupancies must be paired and length ≥ 3")
    if np.any(radii <= 0) or np.any(occupancies <= 0):
        raise ValueError("radii and occupancies must be positive")
    response = np.log(occupancies.astype(np.float64))
    exp_intercept, exp_slope, exp_adj = _fit_affine(radii, response)
    poly_intercept, poly_slope, poly_adj = _fit_affine(np.log(radii), response)
    span = float(radii.max() / radii.min())
    if math.isnan(exp_adj) or math.isnan(poly_adj):
        winner = "undecided"
        margin = math.nan
    elif exp_adj > poly_adj:
        winner = "exponential"
        margin = exp_adj - poly_adj
    elif poly_adj > exp_adj:
        winner = "polynomial"
        margin = poly_adj - exp_adj
    else:
        winner = "undecided"
        margin = 0.0
    return {
        "shells": int(radii.size),
        "span_ratio": span,
        "exponential": {
            "slope": exp_slope,
            "intercept": exp_intercept,
            "adjusted_r_squared": exp_adj,
        },
        "polynomial": {
            "degree": poly_slope,
            "intercept": poly_intercept,
            "adjusted_r_squared": poly_adj,
        },
        "winner": winner,
        "adjusted_r_squared_margin": margin,
    }


def hellinger_affinity(
    r: float,
    d: float,
    n_win: float,
    *,
    steps: int = 256,
) -> float:
    """Hellinger affinity of the endpoint-matched Poisson processes.

    Intensities are normalized so both put mass ``n_win`` on [1, r].
    Returns ∏-equivalent exp(-½ ∫(√λ-√ν)²), i.e. 1 - H².
    """
    if r <= 1.0 or d <= 0.0 or n_win <= 0.0:
        raise ValueError("r>1, d>0, n_win>0 required")
    log_r = math.log(r)
    h = d * log_r / (r - 1.0)  # scaled: ρ_min = 1
    # Unnormalized intensities on t ∈ [1, r].
    # exp: h * exp(h (t-1)); poly: (d/t) * t^d  after matching n(1)=1, n(r)=r^d.
    # Both integrate to r^d - 1. Scale so the integral is n_win.
    t = np.linspace(1.0, r, steps)
    dt = (r - 1.0) / (steps - 1)
    lam_e = h * np.exp(h * (t - 1.0))
    lam_p = d * np.power(t, d - 1.0)
    mass = r**d - 1.0
    scale = n_win / mass
    lam_e = scale * lam_e
    lam_p = scale * lam_p
    integrand = (np.sqrt(lam_e) - np.sqrt(lam_p)) ** 2
    trapz = getattr(np, "trapezoid", None) or np.trapz
    integral = float(trapz(integrand, dx=dt))
    return math.exp(-0.5 * integral)


def lecam_risk_lower_bound(r: float, d: float, n_win: float) -> float:
    """Theorem 2.1: ½(1 - √(2H² - H⁴)) from the Hellinger affinity."""
    affinity = hellinger_affinity(r, d, n_win)
    h2 = max(0.0, min(1.0, 1.0 - affinity))
    tv_upper = math.sqrt(max(0.0, 2.0 * h2 - h2**2))
    return 0.5 * (1.0 - min(1.0, tv_upper))


def measurability(
    *,
    n_win: int,
    span_ratio: float,
    shells: int,
    d_star: float = DEFAULT_D_STAR,
    alpha: float = DEFAULT_ALPHA,
    k_min: int = DEFAULT_K_MIN,
) -> dict[str, object]:
    """Three-way Layer 0 verdict for one occupancy window."""
    if span_ratio <= 1.0 or n_win <= 0 or shells <= 0:
        return {
            "verdict": "unmeasurable",
            "reason": "window is degenerate (need r>1, N>0, k>0)",
            "span_ratio": span_ratio,
            "n_win": n_win,
            "shells": shells,
            "delta": math.nan,
            "lecam_risk_lower_bound": math.nan,
            "critical_n": math.nan,
        }
    delta = span_information(span_ratio, d_star)
    risk = lecam_risk_lower_bound(span_ratio, d_star, float(n_win))
    critical_n = critical_sample_size(span_ratio, d_star)
    if shells < k_min:
        verdict = "unmeasurable"
        reason = f"only {shells} resolvable shells; need k ≥ {k_min}"
    elif risk >= alpha:
        verdict = "unmeasurable"
        reason = (
            f"Le Cam lower bound {risk:.3f} ≥ α={alpha} "
            f"at r={span_ratio:.3f}, N={n_win}, d*={d_star}"
        )
    elif n_win < critical_n:
        verdict = "undecided"
        reason = (
            f"simple-hypothesis pair is not proven unmeasurable, "
            f"but N={n_win} < N*={critical_n:.0f} required for "
            f"Δ√N ≥ {Z_DETECT} at this span"
        )
    else:
        verdict = "measurable"
        reason = (
            f"k={shells} ≥ {k_min}, Le Cam bound {risk:.3f} < α={alpha}, "
            f"and N={n_win} ≥ N*={critical_n:.0f}"
        )
    return {
        "verdict": verdict,
        "reason": reason,
        "span_ratio": span_ratio,
        "n_win": n_win,
        "shells": shells,
        "d_star": d_star,
        "alpha": alpha,
        "k_min": k_min,
        "delta": delta,
        "midpoint_ratio": midpoint_ratio(span_ratio, d_star),
        "lecam_risk_lower_bound": risk,
        "critical_n": critical_n,
    }


def analyze_distances(
    distances: Sequence[float],
    *,
    d_star: float = DEFAULT_D_STAR,
    alpha: float = DEFAULT_ALPHA,
    k_min: int = DEFAULT_K_MIN,
    min_count: int = 3,
    max_fraction: float = 0.5,
    points: int = 32,
    span: float | None = None,
) -> dict[str, object]:
    """Measurability plus, if licensed, a growth-class call."""
    profile = occupancy_profile(
        distances,
        min_count=min_count,
        max_fraction=max_fraction,
        points=points,
    )
    if profile is None:
        return {
            "measurable": False,
            "verdict": "unmeasurable",
            "reason": "insufficient distinct radii for an occupancy profile",
            "growth_class": None,
        }
    radii, occupancies = profile
    if span is not None:
        restricted = restrict_span(radii, occupancies, span)
        if restricted is None:
            return {
                "measurable": False,
                "verdict": "unmeasurable",
                "reason": f"cannot form a window of span {span}",
                "growth_class": None,
            }
        radii, occupancies = restricted
    n_win = int(occupancies[-1] - occupancies[0] + (occupancies[0] > 0))
    span_ratio = float(radii.max() / radii.min())
    gate = classify_occupancy(radii, occupancies)
    status = measurability(
        n_win=max(n_win, int(occupancies[-1])),
        span_ratio=span_ratio,
        shells=int(radii.size),
        d_star=d_star,
        alpha=alpha,
        k_min=k_min,
    )
    growth_class = gate["winner"] if status["verdict"] == "measurable" else None
    return {
        "measurable": status["verdict"] == "measurable",
        "verdict": status["verdict"],
        "reason": status["reason"],
        "growth_class": growth_class,
        "gate": gate,
        "measurability": status,
    }
