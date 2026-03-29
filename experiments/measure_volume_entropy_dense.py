#!/usr/bin/env python3
"""
Dense-window volume entropy measurement.
Uses shorter windows (1.2s) and small hop (0.3s) to get ~300-400 points
on the SPD manifold, enabling proper exponential growth estimation.
"""
import json, sys, time, numpy as np
from numpy.linalg import eigh
from pathlib import Path

LN2 = np.log(2)

def mat_log(C):
    w, V = eigh(C)
    w = np.clip(w, 1e-10, None)
    return V @ np.diag(np.log(w)) @ V.T

def log_euclidean_distance(L1, L2):
    D = L1 - L2
    return float(np.sqrt(np.sum(D * D)))

def distance_matrix_loge(log_covs):
    n = len(log_covs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = log_euclidean_distance(log_covs[i], log_covs[j])
            D[i, j] = D[j, i] = d
    return D

def tri_kappa_bootstrap(D, ns=2000, B=500, seed=0):
    n = D.shape[0]
    if n < 3:
        return float('nan'), (float('nan'), float('nan'))
    med = np.median(D[np.triu_indices(n, 1)])
    if med > 1e-10:
        Dn = D / med
    else:
        Dn = D.copy()
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
    return float(np.median(ests)), (float(np.percentile(ests, 2.5)), float(np.percentile(ests, 97.5)))

def estimate_volume_entropy(D, n_scales=30):
    n_pts = D.shape[0]
    if n_pts < 20:
        return {"h_vol_nats": float('nan'), "error": "too few points"}
    slopes_nats = []
    for center in range(min(n_pts, 50)):
        dists = np.sort(D[center, :])
        dists = dists[dists > 0]
        if len(dists) < 10:
            continue
        unique_d = np.unique(dists)
        counts = np.array([np.sum(dists <= r) for r in unique_d])
        mask = counts > 0
        R = unique_d[mask]
        logN = np.log(counts[mask].astype(float))
        if len(R) < 5:
            continue
        lo, hi = int(0.2 * len(R)), int(0.8 * len(R))
        if hi - lo < 3:
            lo, hi = 0, len(R)
        R_mid = R[lo:hi]
        logN_mid = logN[lo:hi]
        A = np.vstack([R_mid, np.ones_like(R_mid)]).T
        result = np.linalg.lstsq(A, logN_mid, rcond=None)
        slope = result[0][0]
        ss_res = np.sum((logN_mid - (slope * R_mid + result[0][1])) ** 2)
        ss_tot = np.sum((logN_mid - logN_mid.mean()) ** 2)
        r2 = 1 - ss_res / max(ss_tot, 1e-30)
        if r2 > 0.5:
            slopes_nats.append(slope)
    if not slopes_nats:
        return {"h_vol_nats": float('nan'), "error": "no good fits"}
    arr = np.array(slopes_nats)
    return {
        "h_vol_nats": float(np.median(arr)),
        "h_vol_bits": float(np.median(arr) / LN2),
        "h_vol_nats_std": float(np.std(arr)),
        "n_good_fits": len(arr),
    }

def correlation_dimension(D, n_scales=25):
    n = D.shape[0]
    triu_idx = np.triu_indices(n, 1)
    dists = D[triu_idx]
    pos = dists[dists > 0]
    if len(pos) < 20:
        return {"d_corr": float('nan')}
    d_min = np.percentile(pos, 10)
    d_max = np.percentile(pos, 90)
    if d_min >= d_max or d_min <= 0:
        return {"d_corr": float('nan')}
    radii = np.logspace(np.log10(d_min), np.log10(d_max), n_scales)
    C_r = np.array([np.sum(dists < r) / len(dists) for r in radii])
    mask = C_r > 0
    if mask.sum() < 4:
        return {"d_corr": float('nan')}
    log_r = np.log(radii[mask])
    log_C = np.log(C_r[mask])
    A_mat = np.vstack([log_r, np.ones_like(log_r)]).T
    result = np.linalg.lstsq(A_mat, log_C, rcond=None)
    d_corr = result[0][0]
    local = np.diff(log_C) / np.diff(log_r) if len(log_r) > 1 else [d_corr]
    return {"d_corr": round(float(d_corr), 3), "local_median": round(float(np.median(local)), 3)}

def select_stable_neurons(spks, cap, fano_max):
    tot = spks.sum(axis=(1, 2))
    cf = spks.reshape(spks.shape[0], -1)
    means = cf.mean(axis=1)
    vars_ = cf.var(axis=1)
    fano = np.where(means > 0, vars_ / np.maximum(means, 1e-12), np.inf)
    T = cf.shape[1]
    half = T // 2
    has_both = (cf[:, :half].sum(axis=1) > 0) & (cf[:, half:].sum(axis=1) > 0)
    idx = np.where((fano < fano_max) & has_both)[0]
    if len(idx) == 0:
        idx = np.arange(spks.shape[0])
    idx_sorted = idx[np.argsort(tot[idx])]
    return idx_sorted[-min(cap, len(idx_sorted)):]

CONFIGS = [
    {"name": "dense_1.2s",    "window_s": 1.2, "hop_s": 0.3, "bin_s": 0.3},
    {"name": "dense_2.4s",    "window_s": 2.4, "hop_s": 0.6, "bin_s": 0.3},
    {"name": "original_4.8s", "window_s": 4.8, "hop_s": 2.4, "bin_s": 0.3},
]

NEURON_CAP = 180
FANO_MAX = 2.0

def process_session(sp, config):
    t0 = time.time()
    name = sp.stem
    arr = np.load(sp, allow_pickle=True)
    dat = arr["dat"][0]
    spks_full = dat["spks"]
    N_all, Ttr, B = spks_full.shape
    bin_size = 0.01

    keep = select_stable_neurons(spks_full, NEURON_CAP, FANO_MAX)
    spks = spks_full[keep]
    N = spks.shape[0]
    if N < 10:
        return {"session": name, "config": config["name"], "error": "neurons=%d" % N}

    X_flat = spks.reshape(N, -1).T
    rebin = max(1, int(round(config["bin_s"] / bin_size)))
    n_bins = X_flat.shape[0] // rebin
    X = X_flat[:n_bins * rebin].reshape(n_bins, rebin, N).sum(axis=1)

    win = max(3, int(round(config["window_s"] / config["bin_s"])))
    hop = max(1, int(round(config["hop_s"] / config["bin_s"])))
    covs = []
    log_covs = []
    for i in range(0, n_bins - win + 1, hop):
        seg = X[i:i+win]
        if len(seg) >= 3:
            C = np.cov(seg.T) + 1e-6 * np.eye(N)
            covs.append(C)
            log_covs.append(mat_log(C))

    nw = len(log_covs)
    if nw < 20:
        return {"session": name, "config": config["name"], "error": "windows=%d" % nw, "n_neurons": N}

    D = distance_matrix_loge(log_covs)
    k, k_ci = tri_kappa_bootstrap(D)
    vol = estimate_volume_entropy(D)
    cdim = correlation_dimension(D)

    h_pred = float(np.sqrt(k) / LN2) if k > 0 else float('nan')
    h_meas = vol.get("h_vol_bits", float('nan'))
    h_nats = vol.get("h_vol_nats", float('nan'))
    n_impl = 1 + h_nats / np.sqrt(max(k, 1e-10)) if (k > 0 and np.isfinite(h_nats)) else float('nan')

    # Tangent space PCA
    d = log_covs[0].shape[0]
    idx = np.triu_indices(d)
    V = np.array([L[idx] for L in log_covs])
    V = V - V.mean(axis=0)
    _, s, _ = np.linalg.svd(V, full_matrices=False)
    s2 = s ** 2
    n_pr = float((s2.sum())**2 / (s2**2).sum()) if s2.sum() > 0 else 0
    cumvar = np.cumsum(s2) / max(s2.sum(), 1e-30)

    return {
        "session": name, "config": config["name"], "n_neurons": N, "n_windows": nw,
        "kappa": round(float(k), 4),
        "kappa_ci": [round(float(k_ci[0]), 4), round(float(k_ci[1]), 4)],
        "h_vol_bits": round(h_meas, 4) if np.isfinite(h_meas) else None,
        "h_predicted_n2": round(h_pred, 4),
        "n_implied": round(float(n_impl), 3) if np.isfinite(n_impl) else None,
        "d_corr": cdim.get("d_corr"),
        "traj_n_pr": round(n_pr, 2),
        "traj_80pct": int(np.searchsorted(cumvar, 0.80)) + 1,
        "elapsed_s": round(time.time() - t0, 1),
    }

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sessions", type=int, default=5)
    args = parser.parse_args()

    sessions = sorted(Path(args.data_dir).glob("steinmetz_session_*.npz"))[:args.sessions]
    print("Dense Volume Entropy -- %d sessions x %d configs" % (len(sessions), len(CONFIGS)))
    print("Prediction for n=2: h = sqrt(kappa)/ln2")
    print("=" * 70)

    all_results = []
    for i, sp in enumerate(sessions):
        print("\n[%d/%d] %s" % (i+1, len(sessions), sp.name))
        for cfg in CONFIGS:
            r = process_session(sp, cfg)
            cn = cfg["name"]
            if "error" not in r:
                hv = r["h_vol_bits"] if r["h_vol_bits"] is not None else "N/A"
                hp = r["h_predicted_n2"]
                ni = r["n_implied"] if r["n_implied"] is not None else "N/A"
                dc = r["d_corr"] if r["d_corr"] is not None else "N/A"
                tp = r["traj_n_pr"]
                print("  %-15s | wins=%3d | k=%.4f | h_vol=%s | h_pred(n=2)=%.4f | n_impl=%s | d_corr=%s | PR=%.1f" %
                      (cn, r["n_windows"], r["kappa"], hv, hp, ni, dc, tp))
            else:
                print("  %-15s | ERROR: %s" % (cn, r["error"]))
            all_results.append(r)

    # Summary by config
    print("\n" + "=" * 70)
    print("SUMMARY BY WINDOW CONFIGURATION")
    print("=" * 70)
    for cfg in CONFIGS:
        valid = [r for r in all_results if r.get("config") == cfg["name"] and "error" not in r]
        if valid:
            k_mean = np.mean([r["kappa"] for r in valid])
            h_vols = [r["h_vol_bits"] for r in valid if r.get("h_vol_bits") is not None]
            h_preds = [r["h_predicted_n2"] for r in valid]
            n_impls = [r["n_implied"] for r in valid if r.get("n_implied") is not None]
            d_corrs = [r["d_corr"] for r in valid if r.get("d_corr") is not None and np.isfinite(r["d_corr"])]
            nws = [r["n_windows"] for r in valid]
            print("\n%s:" % cfg["name"])
            print("  kappa = %.4f" % k_mean)
            if h_vols:
                print("  h_vol = %.4f +/- %.4f bits" % (np.mean(h_vols), np.std(h_vols)))
            if h_preds:
                print("  h_pred(n=2) = %.4f" % np.mean(h_preds))
            if n_impls:
                print("  n_implied = %.3f +/- %.3f" % (np.mean(n_impls), np.std(n_impls)))
            if d_corrs:
                print("  d_corr = %.2f +/- %.2f" % (np.mean(d_corrs), np.std(d_corrs)))
            print("  n_windows = %.0f" % np.mean(nws))

    with open(args.output, "w") as f:
        json.dump({"configs": CONFIGS, "results": all_results}, f, indent=2, default=str)
    print("\nWrote %s" % args.output)

if __name__ == "__main__":
    main()
