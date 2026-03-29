#!/usr/bin/env python3
"""
Peri-Stimulus κ(t) and h(t) Dynamics.

Tests the field equation prediction: κ̇ = −Aκ² + B(t)κ

By pooling covariance matrices across trials at each peri-stimulus
time point, we track κ(t) and h(t) around stimulus onset. The theory
predicts:
  - κ shows transient perturbation, then exponential relaxation to κ*
  - h peaks during stimulus processing, then decays
  - n(t) = 1 + h·ln2/√κ remains constant at ~2 throughout

Usage:
    python peri_stimulus_dynamics.py \
        --data-dir /path/to/steinmetz_cache/ \
        --output peri_stimulus_results.json \
        [--sessions 5]
"""
import json, time, argparse, numpy as np
from numpy.linalg import eigh
from pathlib import Path
from scipy.optimize import curve_fit

LN2 = np.log(2)

# ── SPD geometry (standalone) ────────────────────────────────────────────

def mat_log(C):
    w, V = eigh(C)
    w = np.clip(w, 1e-10, None)
    return V @ np.diag(np.log(w)) @ V.T

def distance_matrix_loge(log_covs):
    n = len(log_covs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            diff = log_covs[i] - log_covs[j]
            D[i, j] = D[j, i] = float(np.sqrt(np.sum(diff * diff)))
    return D

def tri_kappa_bootstrap(D, ns=2000, B=200, seed=0):
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
    return float(np.median(ests)), (float(np.percentile(ests, 2.5)),
                                     float(np.percentile(ests, 97.5)))

def estimate_volume_entropy(D, n_centers=40, r2_threshold=0.5):
    n_pts = D.shape[0]
    if n_pts < 15:
        return {"h_vol_nats": float('nan'), "error": "too few points"}
    slopes = []
    for center in range(min(n_pts, n_centers)):
        dists = np.sort(D[center, :])
        dists = dists[dists > 0]
        if len(dists) < 8:
            continue
        unique_d = np.unique(dists)
        counts = np.array([np.sum(dists <= r) for r in unique_d])
        mask = counts > 0
        R = unique_d[mask]
        logN = np.log(counts[mask].astype(float))
        if len(R) < 4:
            continue
        lo, hi = int(0.2 * len(R)), int(0.8 * len(R))
        if hi - lo < 3:
            lo, hi = 0, len(R)
        R_mid, logN_mid = R[lo:hi], logN[lo:hi]
        A = np.vstack([R_mid, np.ones_like(R_mid)]).T
        result = np.linalg.lstsq(A, logN_mid, rcond=None)
        slope = result[0][0]
        ss_res = np.sum((logN_mid - (slope * R_mid + result[0][1])) ** 2)
        ss_tot = np.sum((logN_mid - logN_mid.mean()) ** 2)
        r2 = 1 - ss_res / max(ss_tot, 1e-30)
        if r2 > r2_threshold:
            slopes.append(slope)
    if not slopes:
        return {"h_vol_nats": float('nan'), "error": "no good fits"}
    arr = np.array(slopes)
    return {"h_vol_nats": float(np.median(arr)), "h_vol_bits": float(np.median(arr) / LN2),
            "h_vol_std": float(np.std(arr)), "n_fits": len(arr)}

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

# ── Peri-stimulus analysis ───────────────────────────────────────────────

NEURON_CAP = 180
FANO_MAX = 2.0
BIN_MS = 10          # Original Steinmetz bin size
REBIN_MS = 300       # Re-bin to 300ms for covariance stability
WINDOW_MS = 600      # Covariance window: 600ms (2 re-bins)
# Time grid around stimulus onset
T_PRE_MS = -1500     # 1.5s before stimulus
T_POST_MS = 2700     # 2.7s after stimulus
T_STEP_MS = 300      # 300ms steps


def exponential_decay(t, kappa_star, delta_k0, tau):
    """Exponential relaxation model: κ(t) = κ* + δκ₀ · exp(−t/τ)"""
    return kappa_star + delta_k0 * np.exp(-t / max(tau, 1e-6))


def process_session(sp):
    t0 = time.time()
    name = sp.stem
    arr = np.load(sp, allow_pickle=True)
    dat = arr["dat"][0]
    spks_full = dat["spks"]          # (N_neurons, N_trials, N_bins)
    N_all, Ttr, B = spks_full.shape
    bin_s = 0.01                      # 10ms bins

    # Trial metadata
    stim_onset = float(dat.get("stim_onset", 0.5))  # seconds into trial
    contrast_left = dat.get("contrast_left", np.zeros(Ttr))
    contrast_right = dat.get("contrast_right", np.zeros(Ttr))
    feedback_type = dat.get("feedback_type", np.zeros(Ttr))
    mouse = str(dat.get("mouse_name", "unknown"))

    # Select neurons
    keep = select_stable_neurons(spks_full, NEURON_CAP, FANO_MAX)
    spks = spks_full[keep]   # (N, Ttr, B)
    N = spks.shape[0]
    if N < 10:
        return {"session": name, "error": "too few neurons (%d)" % N}

    # Stimulus onset in bin units
    stim_bin = int(stim_onset / bin_s)
    rebin_factor = REBIN_MS // BIN_MS  # 30 bins per re-bin

    # Define time grid (ms relative to stimulus)
    time_grid_ms = list(range(T_PRE_MS, T_POST_MS + 1, T_STEP_MS))
    win_rebins = WINDOW_MS // REBIN_MS  # 2 re-bins per window

    # For each time point: collect one covariance per trial
    dynamics = {"time_ms": time_grid_ms, "kappa": [], "kappa_ci": [],
                "h_vol_bits": [], "n_implied": [], "n_trials_used": []}
    dynamics_high = {"kappa": [], "h_vol_bits": [], "n_implied": []}
    dynamics_low = {"kappa": [], "h_vol_bits": [], "n_implied": []}

    # Classify trials by contrast
    max_contrast = np.maximum(contrast_left, contrast_right)
    high_mask = max_contrast > 0.5
    low_mask = (max_contrast > 0) & (max_contrast <= 0.5)

    for t_ms in time_grid_ms:
        # Center of window in original bins, relative to trial start
        center_bin = stim_bin + int(t_ms / BIN_MS)
        # Extract re-binned window
        start_bin = center_bin - (win_rebins * rebin_factor) // 2
        end_bin = start_bin + win_rebins * rebin_factor

        if start_bin < 0 or end_bin > B:
            dynamics["kappa"].append(None)
            dynamics["kappa_ci"].append(None)
            dynamics["h_vol_bits"].append(None)
            dynamics["n_implied"].append(None)
            dynamics["n_trials_used"].append(0)
            dynamics_high["kappa"].append(None)
            dynamics_high["h_vol_bits"].append(None)
            dynamics_high["n_implied"].append(None)
            dynamics_low["kappa"].append(None)
            dynamics_low["h_vol_bits"].append(None)
            dynamics_low["n_implied"].append(None)
            continue

        # Collect one covariance matrix per trial
        log_covs_all = []
        log_covs_high = []
        log_covs_low = []

        for trial in range(Ttr):
            seg = spks[:, trial, start_bin:end_bin]  # (N, window_bins)
            # Re-bin: sum over rebin_factor bins
            n_rebins = seg.shape[1] // rebin_factor
            if n_rebins < 2:
                continue
            rebinned = seg[:, :n_rebins * rebin_factor].reshape(N, n_rebins, rebin_factor).sum(axis=2)
            # Covariance: (N, N) from (n_rebins, N) transposed
            X = rebinned.T  # (n_rebins, N)
            if X.shape[0] < 2:
                continue
            C = np.cov(X.T) + 1e-6 * np.eye(N)
            L = mat_log(C)
            log_covs_all.append(L)
            if high_mask[trial]:
                log_covs_high.append(L)
            elif low_mask[trial]:
                log_covs_low.append(L)

        # Compute κ and h for all trials at this time point
        def compute_kh(log_covs):
            if len(log_covs) < 15:
                return None, None, None, None
            D = distance_matrix_loge(log_covs)
            k, k_ci = tri_kappa_bootstrap(D, ns=1500, B=200)
            vol = estimate_volume_entropy(D)
            h_bits = vol.get("h_vol_bits", float('nan'))
            h_nats = vol.get("h_vol_nats", float('nan'))
            n_impl = 1 + h_nats / np.sqrt(max(k, 1e-10)) if (k > 0 and np.isfinite(h_nats)) else float('nan')
            return k, k_ci, h_bits, n_impl

        k_all, k_ci_all, h_all, n_all = compute_kh(log_covs_all)
        k_hi, _, h_hi, n_hi = compute_kh(log_covs_high)
        k_lo, _, h_lo, n_lo = compute_kh(log_covs_low)

        dynamics["kappa"].append(round(k_all, 5) if k_all is not None else None)
        dynamics["kappa_ci"].append([round(k_ci_all[0], 5), round(k_ci_all[1], 5)] if k_ci_all is not None else None)
        dynamics["h_vol_bits"].append(round(h_all, 4) if h_all is not None and np.isfinite(h_all) else None)
        dynamics["n_implied"].append(round(n_all, 3) if n_all is not None and np.isfinite(n_all) else None)
        dynamics["n_trials_used"].append(len(log_covs_all))

        dynamics_high["kappa"].append(round(k_hi, 5) if k_hi is not None else None)
        dynamics_high["h_vol_bits"].append(round(h_hi, 4) if h_hi is not None and np.isfinite(h_hi) else None)
        dynamics_high["n_implied"].append(round(n_hi, 3) if n_hi is not None and np.isfinite(n_hi) else None)

        dynamics_low["kappa"].append(round(k_lo, 5) if k_lo is not None else None)
        dynamics_low["h_vol_bits"].append(round(h_lo, 4) if h_lo is not None and np.isfinite(h_lo) else None)
        dynamics_low["n_implied"].append(round(n_lo, 3) if n_lo is not None and np.isfinite(n_lo) else None)

    # Relaxation fit (on post-stimulus κ)
    relax_fit = None
    valid_ks = [(t, k) for t, k in zip(time_grid_ms, dynamics["kappa"])
                if k is not None and t >= 0]
    if len(valid_ks) >= 4:
        ts = np.array([t / 1000.0 for t, _ in valid_ks])  # seconds
        ks = np.array([k for _, k in valid_ks])
        try:
            popt, pcov = curve_fit(exponential_decay, ts, ks,
                                   p0=[np.mean(ks), ks[0] - np.mean(ks), 0.5],
                                   maxfev=5000)
            kappa_star, delta_k0, tau = popt
            ss_res = np.sum((ks - exponential_decay(ts, *popt)) ** 2)
            ss_tot = np.sum((ks - ks.mean()) ** 2)
            r2 = 1 - ss_res / max(ss_tot, 1e-30)
            relax_fit = {
                "kappa_star": round(float(kappa_star), 5),
                "delta_kappa_0": round(float(delta_k0), 5),
                "tau_seconds": round(float(tau), 4),
                "r2": round(float(r2), 4),
            }
        except Exception:
            relax_fit = {"error": "fit failed"}

    # n(t) stability
    valid_ns = [n for n in dynamics["n_implied"] if n is not None and np.isfinite(n)]
    n_stability = {
        "mean": round(float(np.mean(valid_ns)), 3) if valid_ns else None,
        "std": round(float(np.std(valid_ns)), 3) if valid_ns else None,
        "cv": round(float(np.std(valid_ns) / np.mean(valid_ns)), 4) if valid_ns and np.mean(valid_ns) > 0 else None,
    }

    return {
        "session": name,
        "mouse": mouse,
        "n_neurons": N,
        "n_trials": Ttr,
        "n_trials_high_contrast": int(high_mask.sum()),
        "n_trials_low_contrast": int(low_mask.sum()),
        "stim_onset_s": stim_onset,
        "dynamics": dynamics,
        "high_contrast": dynamics_high,
        "low_contrast": dynamics_low,
        "relaxation_fit": relax_fit,
        "n_stability": n_stability,
        "elapsed_s": round(time.time() - t0, 1),
    }


def main():
    parser = argparse.ArgumentParser(description="Peri-stimulus κ(t) and h(t) dynamics")
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sessions", type=int, default=None)
    args = parser.parse_args()

    sessions = sorted(Path(args.data_dir).glob("steinmetz_session_*.npz"))
    if args.sessions:
        sessions = sessions[:args.sessions]

    print("Peri-Stimulus Dynamics: %d sessions" % len(sessions))
    print("Field equation: kappa_dot = -A*kappa^2 + B(t)*kappa")
    print("Time grid: %d to %d ms, step %d ms" % (T_PRE_MS, T_POST_MS, T_STEP_MS))
    print("=" * 70)

    all_results = []
    for i, sp in enumerate(sessions):
        print("\n[%d/%d] %s ..." % (i + 1, len(sessions), sp.name), flush=True)
        try:
            r = process_session(sp)
            if "error" not in r:
                ns = r["n_stability"]
                rf = r["relaxation_fit"]
                print("  %s | %d trials | n_mean=%.3f +/- %.3f" %
                      (r["mouse"], r["n_trials"],
                       ns["mean"] if ns["mean"] else 0,
                       ns["std"] if ns["std"] else 0))
                if rf and "tau_seconds" in rf:
                    print("  Relaxation: tau=%.3fs, delta_k=%.5f, R2=%.3f" %
                          (rf["tau_seconds"], rf["delta_kappa_0"], rf["r2"]))
                # Print κ trajectory
                ks = [k for k in r["dynamics"]["kappa"] if k is not None]
                if ks:
                    print("  kappa(t): min=%.4f max=%.4f range=%.5f" %
                          (min(ks), max(ks), max(ks) - min(ks)))
            else:
                print("  ERROR: %s" % r["error"])
            all_results.append(r)
        except Exception as e:
            print("  EXCEPTION: %s" % e)
            import traceback; traceback.print_exc()
            all_results.append({"session": sp.stem, "error": str(e)})

    # Cohort summary
    valid = [r for r in all_results if "error" not in r]
    if valid:
        print("\n" + "=" * 70)
        print("COHORT SUMMARY: %d sessions" % len(valid))
        print("=" * 70)

        # Average κ(t) across sessions
        time_grid = valid[0]["dynamics"]["time_ms"]
        avg_kappa = []
        avg_h = []
        avg_n = []
        for ti in range(len(time_grid)):
            ks = [r["dynamics"]["kappa"][ti] for r in valid
                  if r["dynamics"]["kappa"][ti] is not None]
            hs = [r["dynamics"]["h_vol_bits"][ti] for r in valid
                  if r["dynamics"]["h_vol_bits"][ti] is not None]
            ns = [r["dynamics"]["n_implied"][ti] for r in valid
                  if r["dynamics"]["n_implied"][ti] is not None]
            avg_kappa.append(round(float(np.mean(ks)), 5) if ks else None)
            avg_h.append(round(float(np.mean(hs)), 4) if hs else None)
            avg_n.append(round(float(np.mean(ns)), 3) if ns else None)

        print("\nGrand-average κ(t):")
        for t, k, h, n in zip(time_grid, avg_kappa, avg_h, avg_n):
            marker = " <-- STIM" if t == 0 else ""
            print("  t=%+5dms: κ=%-8s h=%-8s n=%-8s%s" %
                  (t, k if k else "---", h if h else "---", n if n else "---", marker))

        # n(t) stability across cohort
        all_ns = [n for n in avg_n if n is not None]
        if all_ns:
            print("\nn(t) stability: mean=%.3f, std=%.3f, CV=%.4f" %
                  (np.mean(all_ns), np.std(all_ns),
                   np.std(all_ns) / np.mean(all_ns)))

        # Relaxation fits
        taus = [r["relaxation_fit"]["tau_seconds"] for r in valid
                if r.get("relaxation_fit") and "tau_seconds" in r["relaxation_fit"]]
        if taus:
            print("\nRelaxation timescale tau: %.3f +/- %.3f s (n=%d sessions)" %
                  (np.mean(taus), np.std(taus), len(taus)))

    output = {
        "experiment": "peri_stimulus_dynamics",
        "field_equation": "kappa_dot = -A*kappa^2 + B(t)*kappa",
        "config": {
            "rebin_ms": REBIN_MS, "window_ms": WINDOW_MS,
            "t_pre_ms": T_PRE_MS, "t_post_ms": T_POST_MS, "t_step_ms": T_STEP_MS,
            "neuron_cap": NEURON_CAP, "fano_max": FANO_MAX,
        },
        "grand_average": {
            "time_ms": time_grid if valid else [],
            "kappa": avg_kappa if valid else [],
            "h_vol_bits": avg_h if valid else [],
            "n_implied": avg_n if valid else [],
        },
        "sessions": all_results,
    }
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print("\nWrote %s" % args.output)


if __name__ == "__main__":
    main()
