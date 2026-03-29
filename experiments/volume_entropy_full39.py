#!/usr/bin/env python3
"""
Full 39-session volume entropy at 2.4s windows with brain region metadata.
Runs the volume entropy measurement that showed n_implied ≈ 2.0 in the pilot,
and extracts brain region composition for stratified analysis.
"""
import json, sys, time, numpy as np
from numpy.linalg import eigh
from pathlib import Path
from collections import Counter

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

def estimate_volume_entropy(D):
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

# Classify brain regions by connectivity architecture
HIERARCHICAL_REGIONS = {
    'VISp', 'VISl', 'VISrl', 'VISam', 'VISpm', 'VISal',  # visual cortex
    'MOp', 'MOs',  # motor cortex
    'SSp', 'SSs',  # somatosensory cortex
    'ACA',  # anterior cingulate (hierarchical feedback)
    'RSP',  # retrosplenial
}
RECURRENT_REGIONS = {
    'MD', 'VPM', 'VPL', 'POL', 'LP', 'LGd', 'LD',  # thalamic nuclei
    'PL', 'ILA', 'ORB',  # prefrontal (highly recurrent)
    'CA1', 'CA3', 'DG', 'SUB',  # hippocampal (recurrent loops)
}
SUBCORTICAL_REGIONS = {
    'CP', 'GPe', 'GPi', 'SNr', 'STN',  # basal ganglia
    'MB', 'MRN', 'RN', 'SCm', 'SCs', 'SCsg',  # midbrain/brainstem
    'LS', 'LSc', 'LSr', 'MS',  # septal
    'LH', 'ZI', 'PT',  # hypothalamus/zona incerta
}

NEURON_CAP = 180
FANO_MAX = 2.0
WINDOW_S = 2.4
HOP_S = 0.6
BIN_S = 0.3

def classify_session(brain_areas):
    """Classify session by dominant connectivity architecture."""
    counts = Counter(brain_areas)
    total = sum(counts.values())

    hier_frac = sum(counts.get(r, 0) for r in HIERARCHICAL_REGIONS) / total
    recur_frac = sum(counts.get(r, 0) for r in RECURRENT_REGIONS) / total
    sub_frac = sum(counts.get(r, 0) for r in SUBCORTICAL_REGIONS) / total

    if hier_frac > 0.4:
        dominant = "hierarchical"
    elif recur_frac > 0.4:
        dominant = "recurrent"
    elif sub_frac > 0.4:
        dominant = "subcortical"
    else:
        dominant = "mixed"

    return {
        "hierarchical_frac": round(hier_frac, 3),
        "recurrent_frac": round(recur_frac, 3),
        "subcortical_frac": round(sub_frac, 3),
        "dominant_type": dominant,
        "n_areas": len(counts),
        "top_areas": dict(counts.most_common(5)),
    }

def process_session(sp):
    t0 = time.time()
    name = sp.stem
    arr = np.load(sp, allow_pickle=True)
    dat = arr["dat"][0]
    spks_full = dat["spks"]
    N_all, Ttr, B = spks_full.shape
    bin_size = 0.01

    # Extract metadata
    mouse = str(dat.get("mouse_name", "unknown"))
    date = str(dat.get("date_exp", "unknown"))
    brain_areas = dat.get("brain_area", np.array([]))
    if isinstance(brain_areas, np.ndarray) and len(brain_areas) > 0:
        area_list = [str(a) for a in brain_areas]
    else:
        area_list = []

    keep = select_stable_neurons(spks_full, NEURON_CAP, FANO_MAX)
    spks = spks_full[keep]
    N = spks.shape[0]
    if N < 10:
        return {"session": name, "error": "neurons=%d" % N}

    # Classify by brain region
    if area_list:
        selected_areas = [area_list[i] for i in keep if i < len(area_list)]
        region_info = classify_session(selected_areas)
    else:
        selected_areas = []
        region_info = {"dominant_type": "unknown"}

    # Build count matrix
    X_flat = spks.reshape(N, -1).T
    rebin = max(1, int(round(BIN_S / bin_size)))
    n_bins = X_flat.shape[0] // rebin
    X = X_flat[:n_bins * rebin].reshape(n_bins, rebin, N).sum(axis=1)

    # Windowed covariances
    win = max(3, int(round(WINDOW_S / BIN_S)))
    hop = max(1, int(round(HOP_S / BIN_S)))
    log_covs = []
    for i in range(0, n_bins - win + 1, hop):
        seg = X[i:i+win]
        if len(seg) >= 3:
            C = np.cov(seg.T) + 1e-6 * np.eye(N)
            log_covs.append(mat_log(C))

    nw = len(log_covs)
    if nw < 20:
        return {"session": name, "error": "windows=%d" % nw, "n_neurons": N}

    D = distance_matrix_loge(log_covs)
    k, k_ci = tri_kappa_bootstrap(D)
    vol = estimate_volume_entropy(D)

    h_pred = float(np.sqrt(k) / LN2) if k > 0 else float('nan')
    h_meas = vol.get("h_vol_bits", float('nan'))
    h_nats = vol.get("h_vol_nats", float('nan'))
    n_impl = 1 + h_nats / np.sqrt(max(k, 1e-10)) if (k > 0 and np.isfinite(h_nats)) else float('nan')

    return {
        "session": name,
        "mouse": mouse,
        "date": date,
        "n_neurons_total": int(N_all),
        "n_neurons_selected": int(N),
        "n_trials": int(Ttr),
        "duration_s": round(Ttr * B * bin_size, 1),
        "n_windows": nw,
        "mean_firing_rate": round(float(spks.sum() / N / Ttr), 2),
        "region_info": region_info,
        "kappa": round(float(k), 4),
        "kappa_ci": [round(float(k_ci[0]), 4), round(float(k_ci[1]), 4)],
        "h_vol_bits": round(h_meas, 4) if np.isfinite(h_meas) else None,
        "h_predicted_n2": round(h_pred, 4),
        "n_implied": round(float(n_impl), 3) if np.isfinite(n_impl) else None,
        "elapsed_s": round(time.time() - t0, 1),
    }

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sessions = sorted(Path(args.data_dir).glob("steinmetz_session_*.npz"))
    print("Full 39-session Volume Entropy (2.4s windows)")
    print("State equation: kappa = (h*ln2/(n-1))^2")
    print("Prediction: n = 2 for hierarchical brain regions")
    print("Sessions: %d" % len(sessions))
    print("=" * 80)

    all_results = []
    for i, sp in enumerate(sessions):
        print("\n[%d/%d] %s ..." % (i+1, len(sessions), sp.name), flush=True)
        try:
            r = process_session(sp)
            if "error" not in r:
                ri = r["region_info"]
                ni = r["n_implied"] if r["n_implied"] is not None else "N/A"
                hv = r["h_vol_bits"] if r["h_vol_bits"] is not None else "N/A"
                print("  mouse=%s  type=%s  k=%.4f  h_vol=%s  n_impl=%s  wins=%d  [hier=%.0f%% recur=%.0f%% sub=%.0f%%]" %
                      (r["mouse"], ri["dominant_type"], r["kappa"], hv, ni,
                       r["n_windows"], ri.get("hierarchical_frac",0)*100,
                       ri.get("recurrent_frac",0)*100, ri.get("subcortical_frac",0)*100))
            else:
                print("  ERROR: %s" % r["error"])
            all_results.append(r)
        except Exception as e:
            print("  EXCEPTION: %s" % e)
            all_results.append({"session": sp.stem, "error": str(e)})

    # Analysis
    valid = [r for r in all_results if "error" not in r and r.get("n_implied") is not None]
    print("\n" + "=" * 80)
    print("FULL COHORT RESULTS: %d valid sessions" % len(valid))
    print("=" * 80)

    if valid:
        ks = [r["kappa"] for r in valid]
        nis = [r["n_implied"] for r in valid]
        hvs = [r["h_vol_bits"] for r in valid if r["h_vol_bits"] is not None]

        print("\nOverall:")
        print("  kappa = %.4f +/- %.4f" % (np.mean(ks), np.std(ks)))
        print("  h_vol = %.4f +/- %.4f bits" % (np.mean(hvs), np.std(hvs)))
        print("  n_implied = %.3f +/- %.3f" % (np.mean(nis), np.std(nis)))
        print("  n_implied median = %.3f" % np.median(nis))
        print("  n_implied 95%% CI = [%.3f, %.3f]" % (np.percentile(nis, 2.5), np.percentile(nis, 97.5)))

        # One-sample t-test against n=2
        from scipy import stats
        t_stat, p_value = stats.ttest_1samp(nis, 2.0)
        print("  t-test vs n=2: t=%.3f, p=%.4f" % (t_stat, p_value))

        # Stratify by region type
        print("\nStratified by dominant connectivity:")
        for rtype in ["hierarchical", "recurrent", "subcortical", "mixed", "unknown"]:
            subset = [r for r in valid if r["region_info"].get("dominant_type") == rtype]
            if subset:
                sub_nis = [r["n_implied"] for r in subset]
                print("  %s (n=%d): n_implied = %.3f +/- %.3f  [range: %.3f - %.3f]" %
                      (rtype, len(subset), np.mean(sub_nis), np.std(sub_nis),
                       min(sub_nis), max(sub_nis)))

        # Correlation: recurrent_frac vs n_implied
        recur_fracs = [r["region_info"].get("recurrent_frac", 0) for r in valid]
        hier_fracs = [r["region_info"].get("hierarchical_frac", 0) for r in valid]
        if len(recur_fracs) > 5:
            r_recur, p_recur = stats.spearmanr(recur_fracs, nis)
            r_hier, p_hier = stats.spearmanr(hier_fracs, nis)
            print("\n  Spearman(recurrent_frac, n_implied): rho=%.3f, p=%.4f" % (r_recur, p_recur))
            print("  Spearman(hierarchical_frac, n_implied): rho=%.3f, p=%.4f" % (r_hier, p_hier))

        # Distribution test: is n_implied consistent with n=2?
        print("\nDistribution of n_implied:")
        bins = [1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0]
        hist, _ = np.histogram(nis, bins=bins)
        for j in range(len(bins)-1):
            bar = "#" * hist[j]
            print("  [%.2f-%.2f): %2d %s" % (bins[j], bins[j+1], hist[j], bar))

    output = {
        "experiment": "volume_entropy_full39",
        "config": {"window_s": WINDOW_S, "hop_s": HOP_S, "bin_s": BIN_S,
                   "neuron_cap": NEURON_CAP, "fano_max": FANO_MAX},
        "n_sessions_total": len(all_results),
        "n_sessions_valid": len(valid),
        "results": all_results,
    }
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print("\nWrote %s" % args.output)

if __name__ == "__main__":
    main()
