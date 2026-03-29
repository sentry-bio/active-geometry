#!/usr/bin/env python3
"""
Volume entropy measurement for neural SPD trajectories.

The state equation κ = (h·ln2/(n-1))² comes from Manning's theorem, where h
is the VOLUME ENTROPY — the exponential growth rate of geodesic balls:

    h_vol = lim_{R→∞} log(Vol(B(p,R))) / R

For constant negative curvature -κ in dimension n:
    h_vol = (n-1)·√κ  [nats]

Previous experiments (measure_h_n_neural.py) tested 5 entropy candidates:
von Neumann, spike rate, spike marginal, VAR(1) innovation. None gave n=2.
But none of them measured volume entropy.

This script:
1. Re-computes SPD covariance trajectories (fast, ~30s per session)
2. Computes the FULL distance matrix in Log-Euclidean and AIRM metrics
3. Estimates volume entropy from geodesic ball growth: N(R) ~ exp(h_vol · R)
4. Estimates intrinsic dimensionality of the log-covariance trajectory
5. Checks consistency with n=2 and the state equation

Prediction for n=2: h_vol = √κ / ln(2) ≈ 0.925 bits (given κ ≈ 0.411)

Usage:
    python measure_volume_entropy.py \
        --data-dir /path/to/steinmetz_cache/ \
        --output volume_entropy_results.json \
        [--sessions 5]   # quick test on first 5 sessions
"""
import argparse
import json
import sys
import time
import numpy as np
from numpy.linalg import eigh, eigvalsh, inv, cholesky
from pathlib import Path
from scipy.optimize import curve_fit
from typing import List, Tuple, Dict, Optional

LN2 = np.log(2)

# =============================================================================
# SPD geometry
# =============================================================================

def mat_log(C: np.ndarray) -> np.ndarray:
    """Matrix logarithm via eigendecomposition for SPD matrices."""
    w, V = eigh(C)
    w = np.clip(w, 1e-10, None)
    return V @ np.diag(np.log(w)) @ V.T


def mat_sqrt_inv(C: np.ndarray) -> np.ndarray:
    """C^{-1/2} via eigendecomposition."""
    w, V = eigh(C)
    w = np.clip(w, 1e-10, None)
    return V @ np.diag(1.0 / np.sqrt(w)) @ V.T


def log_euclidean_distance(L1: np.ndarray, L2: np.ndarray) -> float:
    """||log(C1) - log(C2)||_F."""
    D = L1 - L2
    return float(np.sqrt(np.sum(D * D)))


def airm_distance(C1: np.ndarray, C2: np.ndarray) -> float:
    """Affine-Invariant Riemannian Metric distance on SPD(n).
    d(C1,C2) = ||log(C1^{-1/2} C2 C1^{-1/2})||_F
    """
    C1_inv_sqrt = mat_sqrt_inv(C1)
    M = C1_inv_sqrt @ C2 @ C1_inv_sqrt
    w = eigvalsh(M)
    w = np.clip(w, 1e-10, None)
    return float(np.sqrt(np.sum(np.log(w) ** 2)))


def distance_matrix_loge(log_covs: List[np.ndarray]) -> np.ndarray:
    """Pairwise Log-Euclidean distance matrix."""
    n = len(log_covs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = log_euclidean_distance(log_covs[i], log_covs[j])
            D[i, j] = D[j, i] = d
    return D


def distance_matrix_airm(covs: List[np.ndarray]) -> np.ndarray:
    """Pairwise AIRM distance matrix."""
    n = len(covs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = airm_distance(covs[i], covs[j])
            D[i, j] = D[j, i] = d
    return D


# =============================================================================
# Triangle-excess curvature (same as original pipeline)
# =============================================================================

def tri_kappa_bootstrap(D: np.ndarray, ns: int = 1500, B: int = 500,
                        seed: int = 0) -> Tuple[float, Tuple[float, float]]:
    """Bootstrap triangle-excess curvature estimator."""
    n = D.shape[0]
    if n < 3:
        return float('nan'), (float('nan'), float('nan'))
    med = np.median(D[np.triu_indices(n, 1)])
    if med > 1e-10:
        Dn = D / med
    else:
        Dn = D
    ests = []
    for b in range(B):
        rng = np.random.default_rng(seed + b)
        ex = []
        for _ in range(ns):
            i, j, k = rng.choice(n, 3, replace=False)
            a, b_, c = sorted([Dn[i, j], Dn[j, k], Dn[i, k]])
            if a + b_ > c:
                ex.append((a + b_ - c) / 2)
        ests.append(float(np.median(ex)) if ex else float('nan'))
    ests = np.array([e for e in ests if np.isfinite(e)])
    if len(ests) == 0:
        return float('nan'), (float('nan'), float('nan'))
    return float(np.median(ests)), (
        float(np.percentile(ests, 2.5)),
        float(np.percentile(ests, 97.5)),
    )


# =============================================================================
# Volume entropy estimation
# =============================================================================

def estimate_volume_entropy(D: np.ndarray) -> Dict:
    """
    Estimate volume entropy from geodesic ball growth.

    For a manifold with volume entropy h_vol, the number of points
    within radius R of a center point grows as:

        N(R) ~ C · exp(h_vol · R)

    So: log(N(R)) ≈ h_vol · R + const

    We estimate h_vol as the slope of log(N(R)) vs R, averaged
    over multiple center points.

    For hyperbolic space H^n with curvature -κ:
        h_vol = (n-1)·√κ  [nats]
        h_vol_bits = h_vol / ln(2) = (n-1)·√κ / ln(2)  [bits]
    """
    n_pts = D.shape[0]
    if n_pts < 6:
        return {'h_vol_nats': float('nan'), 'error': 'too few points'}

    # Use each point as a center, compute cumulative count vs radius
    all_slopes_nats = []
    all_slopes_bits = []
    all_r2 = []

    for center in range(n_pts):
        dists = np.sort(D[center, :])
        dists = dists[dists > 0]  # exclude self-distance

        if len(dists) < 4:
            continue

        # Cumulative count at each unique distance
        unique_dists = np.unique(dists)
        counts = np.array([np.sum(dists <= r) for r in unique_dists])

        # Only use points where count > 0
        mask = counts > 0
        R = unique_dists[mask]
        logN = np.log(counts[mask].astype(float))

        if len(R) < 3:
            continue

        # Linear regression: logN = h_vol * R + const
        A = np.vstack([R, np.ones_like(R)]).T
        result = np.linalg.lstsq(A, logN, rcond=None)
        slope, intercept = result[0]

        # R² for quality assessment
        ss_res = np.sum((logN - (slope * R + intercept)) ** 2)
        ss_tot = np.sum((logN - logN.mean()) ** 2)
        r2 = 1 - ss_res / max(ss_tot, 1e-30)

        all_slopes_nats.append(slope)
        all_slopes_bits.append(slope / LN2)
        all_r2.append(r2)

    if not all_slopes_nats:
        return {'h_vol_nats': float('nan'), 'error': 'no valid centers'}

    slopes_nats = np.array(all_slopes_nats)
    slopes_bits = np.array(all_slopes_bits)
    r2_arr = np.array(all_r2)

    return {
        'h_vol_nats': float(np.median(slopes_nats)),
        'h_vol_nats_mean': float(np.mean(slopes_nats)),
        'h_vol_nats_std': float(np.std(slopes_nats)),
        'h_vol_bits': float(np.median(slopes_bits)),
        'h_vol_bits_mean': float(np.mean(slopes_bits)),
        'h_vol_bits_std': float(np.std(slopes_bits)),
        'r2_median': float(np.median(r2_arr)),
        'r2_mean': float(np.mean(r2_arr)),
        'n_centers_used': len(slopes_nats),
        'n_points': n_pts,
    }


# =============================================================================
# Intrinsic dimensionality of log-covariance trajectory
# =============================================================================

def trajectory_intrinsic_dim(log_covs: List[np.ndarray]) -> Dict:
    """
    Intrinsic dimensionality of the log-covariance trajectory in tangent space.

    Unlike measure_n_neural_state which operates on spike counts, this operates
    directly on the space where kappa is measured — the tangent space of SPD(n).

    If the trajectory lives on a 2D submanifold, then n=2 in the state equation
    is the right value.
    """
    T = len(log_covs)
    d = log_covs[0].shape[0]

    # Vectorize upper-triangular elements (the tangent space representation)
    idx = np.triu_indices(d)
    V = np.array([L[idx] for L in log_covs])  # (T, p)
    V = V - V.mean(axis=0)  # center

    # SVD for eigenvalues
    _, s, _ = np.linalg.svd(V, full_matrices=False)
    s2 = s ** 2

    if s2.sum() < 1e-30:
        return {'n_eff': float('nan'), 'error': 'degenerate'}

    # Participation ratio: n_PR = (Σλ)² / Σλ²
    n_pr = float((s2.sum()) ** 2 / (s2 ** 2).sum())

    # Cumulative variance
    cumvar = np.cumsum(s2) / s2.sum()

    # How many dimensions for X% variance
    n_80 = int(np.searchsorted(cumvar, 0.80)) + 1
    n_90 = int(np.searchsorted(cumvar, 0.90)) + 1
    n_95 = int(np.searchsorted(cumvar, 0.95)) + 1

    # Spectral gap: biggest ratio between consecutive eigenvalues
    if len(s2) > 1:
        nonzero = s2[s2 > 1e-12]
        if len(nonzero) > 1:
            ratios = nonzero[:-1] / nonzero[1:]
            n_gap = int(np.argmax(ratios)) + 1
            gap_ratio = float(ratios.max())
        else:
            n_gap = 1
            gap_ratio = float('inf')
    else:
        n_gap = 1
        gap_ratio = float('inf')

    # Correlation dimension (Grassberger-Procaccia) — from distance matrix
    # This is a more robust intrinsic dimension estimator

    return {
        'n_pr': round(n_pr, 3),
        'n_80pct': n_80,
        'n_90pct': n_90,
        'n_95pct': n_95,
        'n_spectral_gap': n_gap,
        'spectral_gap_ratio': round(gap_ratio, 2),
        'n_points': T,
        'ambient_dim': V.shape[1],
        'top_eigenvalues': s2[:10].tolist(),
        'cumulative_variance': cumvar[:min(10, len(cumvar))].tolist(),
    }


def correlation_dimension(D: np.ndarray, n_scales: int = 20) -> Dict:
    """
    Correlation dimension from distance matrix (Grassberger-Procaccia).

    C(r) = (2 / N(N-1)) × #{pairs with d(i,j) < r}
    d_corr = d(log C(r)) / d(log r) in the scaling regime

    This is more robust than PCA for estimating intrinsic dimension
    of a curved manifold.
    """
    n = D.shape[0]
    if n < 6:
        return {'d_corr': float('nan'), 'error': 'too few points'}

    # Extract upper triangle distances
    triu_idx = np.triu_indices(n, 1)
    dists = D[triu_idx]
    n_pairs = len(dists)

    if n_pairs < 10:
        return {'d_corr': float('nan'), 'error': 'too few pairs'}

    # Logarithmic scale of radii
    d_min = np.percentile(dists[dists > 0], 5)
    d_max = np.percentile(dists, 95)
    if d_min >= d_max or d_min <= 0:
        return {'d_corr': float('nan'), 'error': 'degenerate distances'}

    radii = np.logspace(np.log10(d_min), np.log10(d_max), n_scales)

    # Compute correlation integral
    C_r = np.array([np.sum(dists < r) / n_pairs for r in radii])

    # Only use points where C_r > 0
    mask = C_r > 0
    if mask.sum() < 3:
        return {'d_corr': float('nan'), 'error': 'insufficient scaling range'}

    log_r = np.log(radii[mask])
    log_C = np.log(C_r[mask])

    # Linear regression for slope = correlation dimension
    A = np.vstack([log_r, np.ones_like(log_r)]).T
    result = np.linalg.lstsq(A, log_C, rcond=None)
    d_corr = result[0][0]

    # R² for quality
    ss_res = np.sum((log_C - (d_corr * log_r + result[0][1])) ** 2)
    ss_tot = np.sum((log_C - log_C.mean()) ** 2)
    r2 = 1 - ss_res / max(ss_tot, 1e-30)

    # Local slopes for scaling regime analysis
    if len(log_r) > 2:
        local_slopes = np.diff(log_C) / np.diff(log_r)
    else:
        local_slopes = np.array([d_corr])

    return {
        'd_corr': round(float(d_corr), 3),
        'r2': round(float(r2), 4),
        'local_slopes': local_slopes.tolist(),
        'local_slopes_median': round(float(np.median(local_slopes)), 3),
        'scaling_range': [float(radii[mask][0]), float(radii[mask][-1])],
        'n_scale_points': int(mask.sum()),
    }


# =============================================================================
# Neuron selection and data preparation (from original pipeline)
# =============================================================================

def select_stable_neurons(spks, cap, fano_max):
    tot_spikes = spks.sum(axis=(1, 2))
    counts_flat = spks.reshape(spks.shape[0], -1)
    means = counts_flat.mean(axis=1)
    vars_ = counts_flat.var(axis=1)
    fano = np.where(means > 0, vars_ / np.maximum(means, 1e-12), np.inf)
    T = counts_flat.shape[1]
    half = T // 2
    has_both = (
        (counts_flat[:, :half].sum(axis=1) > 0)
        & (counts_flat[:, half:].sum(axis=1) > 0)
    )
    idx = np.where((fano < fano_max) & has_both)[0]
    if len(idx) == 0:
        idx = np.arange(spks.shape[0])
    idx_sorted = idx[np.argsort(tot_spikes[idx])]
    return idx_sorted[-min(cap, len(idx_sorted)):]


# =============================================================================
# Main pipeline
# =============================================================================

BIN_S = 0.3
WINDOW_S = 4.8
HOP_FRACTION = 0.5
NEURON_CAP = 180
FANO_MAX = 2.0
MAX_WINDOWS = 120
NSAMPLES = 1500
BOOTSTRAP = 500


def process_session(session_path: Path) -> Dict:
    """Process one Steinmetz session: volume entropy + intrinsic dimension."""
    t0 = time.time()
    session_name = session_path.stem

    # Load
    arr = np.load(session_path, allow_pickle=True)
    dat = arr['dat'][0]
    spks_full = dat['spks']
    N_all, Ttr, B = spks_full.shape
    bin_size = 0.01  # 10ms bins

    # Select neurons
    keep = select_stable_neurons(spks_full, cap=NEURON_CAP, fano_max=FANO_MAX)
    spks = spks_full[keep]
    N = spks.shape[0]

    if N < 10:
        return {'session': session_name, 'error': f'too few neurons ({N})'}

    # Build count matrix at BIN_S resolution
    # Flatten trials into continuous time
    X_flat = spks.reshape(N, -1).T  # (total_bins, N)
    total_time = Ttr * B * bin_size

    # Re-bin at BIN_S
    orig_bins = X_flat.shape[0]
    rebin_factor = max(1, int(round(BIN_S / bin_size)))
    n_new_bins = orig_bins // rebin_factor
    X_rebinned = X_flat[:n_new_bins * rebin_factor].reshape(n_new_bins, rebin_factor, N).sum(axis=1)

    W = min(MAX_WINDOWS, n_new_bins)
    X_full = X_rebinned[:W]

    # Windowed covariances
    win_bins = max(3, int(round(WINDOW_S / BIN_S)))
    hop = max(1, int(win_bins * HOP_FRACTION))
    covs = []
    log_covs = []
    for i in range(0, W - win_bins + 1, hop):
        seg = X_full[i:i + win_bins]
        if len(seg) > 2:
            C = np.cov(seg.T) + 1e-6 * np.eye(N)
            covs.append(C)
            log_covs.append(mat_log(C))

    n_windows = len(log_covs)
    if n_windows < 6:
        return {'session': session_name, 'error': f'too few windows ({n_windows})',
                'n_neurons': N}

    # ---- Distance matrices ----
    D_loge = distance_matrix_loge(log_covs)

    # AIRM is much slower but geometrically correct
    # Only compute for moderate window counts
    D_airm = None
    if n_windows <= 60:
        try:
            D_airm = distance_matrix_airm(covs)
        except Exception as e:
            D_airm = None

    # ---- Kappa (triangle excess) ----
    k_loge, k_loge_ci = tri_kappa_bootstrap(D_loge, ns=NSAMPLES, B=BOOTSTRAP)

    k_airm = float('nan')
    k_airm_ci = (float('nan'), float('nan'))
    if D_airm is not None:
        k_airm, k_airm_ci = tri_kappa_bootstrap(D_airm, ns=NSAMPLES, B=BOOTSTRAP)

    # ---- Volume entropy (the NEW measurement) ----
    vol_ent_loge = estimate_volume_entropy(D_loge)
    vol_ent_airm = estimate_volume_entropy(D_airm) if D_airm is not None else None

    # ---- Intrinsic dimensionality of log-covariance trajectory ----
    traj_dim = trajectory_intrinsic_dim(log_covs)

    # ---- Correlation dimension ----
    corr_dim_loge = correlation_dimension(D_loge)
    corr_dim_airm = correlation_dimension(D_airm) if D_airm is not None else None

    # ---- State equation check with n=2 ----
    # Prediction: h_vol = √κ [nats] for n=2
    # Or equivalently: h_vol_bits = √κ / ln(2)
    h_predicted_n2_bits = float(np.sqrt(k_loge) / LN2) if k_loge > 0 else float('nan')
    h_predicted_n2_nats = float(np.sqrt(k_loge)) if k_loge > 0 else float('nan')

    # Does the measured h_vol match the n=2 prediction?
    h_vol_measured_nats = vol_ent_loge.get('h_vol_nats', float('nan'))
    h_vol_measured_bits = vol_ent_loge.get('h_vol_bits', float('nan'))

    if np.isfinite(h_vol_measured_nats) and np.isfinite(h_predicted_n2_nats):
        h_ratio = h_vol_measured_nats / h_predicted_n2_nats
        # If ratio ≈ 1: n=2 is correct
        # If ratio > 1: need n > 2
        # n_implied from volume entropy: h_vol = (n-1)·√κ → n = 1 + h_vol/√κ
        n_implied_vol = 1 + h_vol_measured_nats / np.sqrt(max(k_loge, 1e-10))
    else:
        h_ratio = float('nan')
        n_implied_vol = float('nan')

    # Also check AIRM
    if D_airm is not None and vol_ent_airm is not None:
        h_vol_airm_nats = vol_ent_airm.get('h_vol_nats', float('nan'))
        h_predicted_airm_n2 = float(np.sqrt(k_airm)) if k_airm > 0 else float('nan')
        if np.isfinite(h_vol_airm_nats) and np.isfinite(h_predicted_airm_n2):
            n_implied_airm = 1 + h_vol_airm_nats / np.sqrt(max(k_airm, 1e-10))
        else:
            n_implied_airm = float('nan')
    else:
        n_implied_airm = float('nan')

    result = {
        'session': session_name,
        'n_neurons': N,
        'n_windows': n_windows,

        # Curvature
        'kappa_loge': float(k_loge),
        'kappa_loge_ci': [float(k_loge_ci[0]), float(k_loge_ci[1])],
        'kappa_airm': float(k_airm),
        'kappa_airm_ci': [float(k_airm_ci[0]), float(k_airm_ci[1])],

        # Volume entropy (THE KEY MEASUREMENT)
        'volume_entropy_loge': vol_ent_loge,
        'volume_entropy_airm': vol_ent_airm,

        # State equation check (n=2)
        'n2_check': {
            'h_predicted_n2_bits': h_predicted_n2_bits,
            'h_predicted_n2_nats': h_predicted_n2_nats,
            'h_vol_measured_bits': h_vol_measured_bits,
            'h_vol_measured_nats': h_vol_measured_nats,
            'h_ratio': float(h_ratio),  # should be ≈1 for n=2
            'n_implied_from_vol_entropy': float(n_implied_vol),
            'n_implied_from_vol_entropy_airm': float(n_implied_airm),
        },

        # Intrinsic dimension
        'trajectory_dimension': traj_dim,
        'correlation_dimension_loge': corr_dim_loge,
        'correlation_dimension_airm': corr_dim_airm,

        'elapsed_s': round(time.time() - t0, 1),
    }

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Volume entropy measurement for neural SPD trajectories')
    parser.add_argument('--data-dir', required=True,
                        help='Path to steinmetz_cache/ with NPZ files')
    parser.add_argument('--output', required=True,
                        help='Output JSON path')
    parser.add_argument('--sessions', type=int, default=None,
                        help='Max sessions (default: all)')
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    sessions = sorted(data_dir.glob('steinmetz_session_*.npz'))
    if not sessions:
        print(f"No steinmetz_session_*.npz in {data_dir}")
        sys.exit(1)

    if args.sessions:
        sessions = sessions[:args.sessions]

    print(f"Volume Entropy Measurement — {len(sessions)} sessions")
    print(f"Frozen params: bin_s={BIN_S}, window_s={WINDOW_S}, cap={NEURON_CAP}")
    print(f"State equation: κ = (h·ln2/(n-1))²")
    print(f"Prediction for n=2: h_vol ≈ 0.925 bits (given κ ≈ 0.411)")
    print("=" * 70)

    all_results = []
    for i, sp in enumerate(sessions):
        print(f"\n[{i+1}/{len(sessions)}] {sp.name} ...", flush=True)
        try:
            result = process_session(sp)
            if 'error' not in result:
                check = result['n2_check']
                traj = result['trajectory_dimension']
                corr = result.get('correlation_dimension_loge', {})
                print(f"  κ_loge={result['kappa_loge']:.4f}  "
                      f"κ_airm={result['kappa_airm']:.4f}")
                print(f"  h_vol(loge)={check['h_vol_measured_bits']:.4f} bits  "
                      f"h_pred(n=2)={check['h_predicted_n2_bits']:.4f} bits  "
                      f"ratio={check['h_ratio']:.3f}")
                print(f"  n_implied(vol)={check['n_implied_from_vol_entropy']:.2f}  "
                      f"n_implied(airm)={check.get('n_implied_from_vol_entropy_airm', '?')}")
                print(f"  traj_dim: PR={traj['n_pr']}  gap={traj['n_spectral_gap']}  "
                      f"80%={traj['n_80pct']}  90%={traj['n_90pct']}")
                print(f"  corr_dim(loge)={corr.get('d_corr', '?')}")
            else:
                print(f"  ERROR: {result['error']}")
            all_results.append(result)
        except Exception as e:
            print(f"  EXCEPTION: {e}")
            import traceback; traceback.print_exc()
            all_results.append({'session': sp.stem, 'error': str(e)})

    # ---- Cohort summary ----
    valid = [r for r in all_results if 'error' not in r]
    if valid:
        summary = {
            'n_sessions': len(valid),
            'kappa_loge': {
                'mean': float(np.mean([r['kappa_loge'] for r in valid])),
                'std': float(np.std([r['kappa_loge'] for r in valid])),
            },
            'kappa_airm': {
                'mean': float(np.mean([r['kappa_airm'] for r in valid
                                       if np.isfinite(r['kappa_airm'])])),
                'std': float(np.std([r['kappa_airm'] for r in valid
                                     if np.isfinite(r['kappa_airm'])])),
            },
            'volume_entropy_bits': {
                'mean': float(np.mean([r['n2_check']['h_vol_measured_bits']
                                       for r in valid
                                       if np.isfinite(r['n2_check']['h_vol_measured_bits'])])),
                'std': float(np.std([r['n2_check']['h_vol_measured_bits']
                                     for r in valid
                                     if np.isfinite(r['n2_check']['h_vol_measured_bits'])])),
            },
            'h_predicted_n2_bits': {
                'mean': float(np.mean([r['n2_check']['h_predicted_n2_bits']
                                       for r in valid
                                       if np.isfinite(r['n2_check']['h_predicted_n2_bits'])])),
            },
            'n_implied_from_volume_entropy': {
                'mean': float(np.mean([r['n2_check']['n_implied_from_vol_entropy']
                                       for r in valid
                                       if np.isfinite(r['n2_check']['n_implied_from_vol_entropy'])])),
                'std': float(np.std([r['n2_check']['n_implied_from_vol_entropy']
                                     for r in valid
                                     if np.isfinite(r['n2_check']['n_implied_from_vol_entropy'])])),
            },
            'trajectory_dimension_pr': {
                'mean': float(np.mean([r['trajectory_dimension']['n_pr']
                                       for r in valid])),
                'std': float(np.std([r['trajectory_dimension']['n_pr']
                                     for r in valid])),
            },
            'correlation_dimension': {
                'mean': float(np.mean([r['correlation_dimension_loge']['d_corr']
                                       for r in valid
                                       if np.isfinite(r['correlation_dimension_loge'].get('d_corr', float('nan')))])),
                'std': float(np.std([r['correlation_dimension_loge']['d_corr']
                                     for r in valid
                                     if np.isfinite(r['correlation_dimension_loge'].get('d_corr', float('nan')))])),
            },
        }

        print("\n" + "=" * 70)
        print("COHORT SUMMARY")
        print("=" * 70)
        print(f"Sessions: {summary['n_sessions']}")
        print(f"\nCurvature:")
        print(f"  κ (LogE): {summary['kappa_loge']['mean']:.4f} ± {summary['kappa_loge']['std']:.4f}")
        print(f"  κ (AIRM): {summary['kappa_airm']['mean']:.4f} ± {summary['kappa_airm']['std']:.4f}")
        print(f"\nVolume Entropy (THE KEY TEST):")
        print(f"  h_vol measured:  {summary['volume_entropy_bits']['mean']:.4f} ± {summary['volume_entropy_bits']['std']:.4f} bits")
        print(f"  h_vol predicted (n=2): {summary['h_predicted_n2_bits']['mean']:.4f} bits")
        print(f"  n implied:       {summary['n_implied_from_volume_entropy']['mean']:.2f} ± {summary['n_implied_from_volume_entropy']['std']:.2f}")
        print(f"\nTrajectory Dimensionality:")
        print(f"  Participation ratio: {summary['trajectory_dimension_pr']['mean']:.2f} ± {summary['trajectory_dimension_pr']['std']:.2f}")
        print(f"  Correlation dimension: {summary['correlation_dimension']['mean']:.2f} ± {summary['correlation_dimension']['std']:.2f}")
        print()

        n_imp = summary['n_implied_from_volume_entropy']['mean']
        if 1.5 < n_imp < 2.5:
            print(">>> RESULT: Volume entropy CONSISTENT with n=2! <<<")
        elif 2.5 <= n_imp < 3.5:
            print(f">>> RESULT: Volume entropy suggests n≈3 (measured {n_imp:.2f})")
        else:
            print(f">>> RESULT: n_implied = {n_imp:.2f} — investigate further")
    else:
        summary = {'error': 'no valid sessions'}

    output = {
        'experiment': 'volume_entropy_neural',
        'state_equation': 'kappa = (h * ln(2) / (n - 1))^2',
        'manning_formula': 'h_vol = (n-1) * sqrt(kappa) [nats]',
        'prediction_n2': 'h_vol ≈ 0.641 nats ≈ 0.925 bits',
        'frozen_params': {
            'bin_s': BIN_S, 'window_s': WINDOW_S, 'neuron_cap': NEURON_CAP,
            'fano_max': FANO_MAX, 'nsamples': NSAMPLES, 'bootstrap': BOOTSTRAP,
        },
        'summary': summary,
        'sessions': all_results,
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote results to {out_path}")


if __name__ == '__main__':
    main()
