#!/usr/bin/env python3
"""
telescope_selfconsistency.py
============================

Self-consistency telescope: run the same GTDB patristic telescope on two
models (compact2 at κ_train=1.25 and phase1b_redux at κ_train=1.39) and
compare where the Pearson-log correlation peaks.

Key question: does the telescope peak track training κ (model artifact) or
stay at a data-intrinsic value (genuine curvature signal)?

  Hypothesis A (data-intrinsic):  both models peak at κ ≈ 1.28–1.35
  Hypothesis B (training-biased): compact2 peaks at ~1.25, redux peaks at ~1.44

Method
------
1. Encode 250 GTDB genomes (telescope manifest) with each model using
   encode_angular_only (c-independent directions).
2. Dense fixed-κ scan (0.70 → 1.70, step 0.02) using within-tree pairs only
   (bac120 genomes vs bac120, ar53 vs ar53) for maximum signal.
3. Compute Spearman ρ and Pearson-log at each κ for both models.
4. Record peak κ and report side-by-side.

Usage
-----
python3 telescope_selfconsistency.py \
    --checkpoint-a  /fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt \
    --checkpoint-b  /home/rohit/phase1b_redux_kappa139/best.pt \
    --telescope-dir /home/rohit/telescope/ \
    --output        /home/rohit/telescope_selfconsistency.json
"""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.stats import pearsonr, spearmanr


# ---------------------------------------------------------------------------
# Poincaré distance
# ---------------------------------------------------------------------------

def poincare_dist_mat(Z, c, eps=1e-7):
    B = Z.shape[0]
    u = Z.unsqueeze(1).expand(B, B, -1)
    v = Z.unsqueeze(0).expand(B, B, -1)
    diff_sq = ((u - v) ** 2).sum(-1)
    u_sq = (u ** 2).sum(-1)
    v_sq = (v ** 2).sum(-1)
    denom = ((1 - c * u_sq) * (1 - c * v_sq)).clamp(min=eps)
    arg = (1 + 2 * c * diff_sq / denom).clamp(min=1 + eps)
    return torch.acosh(arg) / torch.sqrt(c + eps)


# ---------------------------------------------------------------------------
# Load V15Model from checkpoint
# ---------------------------------------------------------------------------

def load_v15model(checkpoint_path, device):
    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = ckpt.get("model_state_dict", ckpt)

    latent_dim  = state["ode_flow.field.0.weight_orig"].shape[1] - 1
    ode_hidden  = state["ode_flow.field.0.bias"].shape[0]
    counts = {
        "bact_fam": state["bact_fam.prototypes"].shape[0],
        "arch_fam": state["arch_fam.prototypes"].shape[0],
        "euk_fam":  state["euk_fam.prototypes"].shape[0],
        "bact_gen": state["bact_gen.prototypes"].shape[0],
        "arch_gen": state["arch_gen.prototypes"].shape[0],
        "euk_gen":  state["euk_gen.prototypes"].shape[0],
    }
    vocab_size = state["encoder.encoder.embed.weight"].shape[0]

    sys.path.insert(0, str(Path.home()))
    sys.path.insert(0, "/fast/sentrybio/scripts")
    from model_v15_5 import V15Model

    model = V15Model(
        vocab_size=vocab_size,
        latent_dim=latent_dim,
        counts=counts,
        ode_hidden=ode_hidden,
    ).to(device)

    state_filtered = {k: v for k, v in state.items() if "curvature_history" not in k}
    model.load_state_dict(state_filtered, strict=False)
    model.eval()
    return model, float(model.live_kappa)


# ---------------------------------------------------------------------------
# Extract embeddings via encode_angular_only
# ---------------------------------------------------------------------------

def extract_embeddings(model, device, manifest_path, max_len=8192, seed=42):
    import random
    random.seed(seed)

    rows = []
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    embeddings = []
    valid_idx  = []
    trees      = []

    with torch.no_grad():
        for i, row in enumerate(rows):
            tp = row.get("tokenized_path", "").strip()
            if not tp or not os.path.exists(tp):
                continue
            try:
                tokens = np.load(tp).astype(np.int64)
                if len(tokens) > max_len:
                    s = random.randint(0, len(tokens) - max_len)
                    tokens = tokens[s : s + max_len]
                t = torch.from_numpy(tokens).unsqueeze(0).to(device)
                z = model.encode_angular_only(t)
                embeddings.append(z.squeeze(0).cpu())
                valid_idx.append(i)
                trees.append(row.get("tree", "unknown"))
            except Exception as e:
                print(f"  Skip {row.get('accession','?')}: {e}")

    return torch.stack(embeddings), valid_idx, trees


# ---------------------------------------------------------------------------
# Dense fixed-κ scan with within-tree mask
# ---------------------------------------------------------------------------

def dense_scan(embeddings, D_pat, trees, kappa_grid):
    """
    For each κ in kappa_grid, compute Spearman ρ and Pearson-log on:
      (a) all pairs
      (b) within-tree pairs (bac120×bac120, ar53×ar53)

    Returns list of dicts.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Z = embeddings.to(device)
    n = Z.shape[0]

    # Build masks
    triu_idx = torch.triu_indices(n, n, offset=1)
    all_mask_np = np.ones(triu_idx.shape[1], dtype=bool)

    within_mask_np = np.array([
        trees[i] == trees[j]
        for i, j in zip(triu_idx[0].numpy(), triu_idx[1].numpy())
    ])

    D_all    = D_pat[triu_idx[0].numpy(), triu_idx[1].numpy()]
    D_within = D_all[within_mask_np]

    results = []
    for kappa in kappa_grid:
        c_t = torch.tensor(kappa, dtype=torch.float32, device=device)
        with torch.no_grad():
            D_hat = poincare_dist_mat(Z, c_t)
        d_flat = D_hat[triu_idx[0], triu_idx[1]].cpu().numpy()

        # All pairs
        rho_all, _ = spearmanr(d_flat, D_all)
        pr_all,  _ = pearsonr(np.log(d_flat + 1e-8), np.log(D_all + 1e-8))

        # Within-tree pairs
        d_within = d_flat[within_mask_np]
        if d_within.shape[0] > 10:
            rho_wt, p_wt = spearmanr(d_within, D_within)
            pr_wt,  _    = pearsonr(np.log(d_within + 1e-8), np.log(D_within + 1e-8))
        else:
            rho_wt = p_wt = pr_wt = float("nan")

        results.append({
            "kappa":          round(float(kappa), 4),
            "spearman_all":   round(float(rho_all), 6),
            "pearson_log_all": round(float(pr_all), 6),
            "spearman_within": round(float(rho_wt), 6),
            "pearson_log_within": round(float(pr_wt), 6),
            "n_within_pairs": int(d_within.shape[0]),
        })
        print(f"  κ={kappa:.3f}  ρ_all={rho_all:.4f}  ρ_within={rho_wt:.4f}"
              f"  pl_all={pr_all:.4f}  pl_within={pr_wt:.4f}")

    return results


def peak_summary(scan_results, label):
    kappas    = [r["kappa"] for r in scan_results]
    rho_all   = [r["spearman_all"] for r in scan_results]
    pl_all    = [r["pearson_log_all"] for r in scan_results]
    rho_wt    = [r["spearman_within"] for r in scan_results]
    pl_wt     = [r["pearson_log_within"] for r in scan_results]

    def safe_peak(vals):
        clean = [(v, k) for v, k in zip(vals, kappas) if math.isfinite(v)]
        if not clean:
            return float("nan"), float("nan")
        best = max(clean, key=lambda x: x[0])
        return best[0], best[1]

    rho_all_val,  rho_all_k   = safe_peak(rho_all)
    pl_all_val,   pl_all_k    = safe_peak(pl_all)
    rho_wt_val,   rho_wt_k    = safe_peak(rho_wt)
    pl_wt_val,    pl_wt_k     = safe_peak(pl_wt)

    print(f"\n  ── {label} peak summary ──")
    print(f"  Spearman(all):     peak ρ={rho_all_val:.4f} at κ={rho_all_k:.4f}")
    print(f"  Pearson-log(all):  peak r={pl_all_val:.4f} at κ={pl_all_k:.4f}")
    print(f"  Spearman(within):  peak ρ={rho_wt_val:.4f} at κ={rho_wt_k:.4f}")
    print(f"  Pearson-log(within): peak r={pl_wt_val:.4f} at κ={pl_wt_k:.4f}")

    return {
        "label": label,
        "spearman_all_peak":    {"rho": rho_all_val, "kappa": rho_all_k},
        "pearson_log_all_peak": {"r":   pl_all_val,  "kappa": pl_all_k},
        "spearman_within_peak": {"rho": rho_wt_val,  "kappa": rho_wt_k},
        "pearson_log_within_peak": {"r": pl_wt_val,  "kappa": pl_wt_k},
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Self-consistency telescope: compact2 vs κ=1.39")
    p.add_argument("--checkpoint-a", default="/fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt",
                   help="Path to model A checkpoint (compact2, κ_train=1.25)")
    p.add_argument("--checkpoint-b", default="/home/rohit/phase1b_redux_kappa139/best.pt",
                   help="Path to model B checkpoint (phase1b_redux, κ_train=1.39)")
    p.add_argument("--telescope-dir", default="/home/rohit/telescope/",
                   help="Directory with telescope_manifest.csv and telescope_patristic.npy")
    p.add_argument("--output", default="/home/rohit/telescope_selfconsistency.json")
    p.add_argument("--max-len", type=int, default=8192)
    p.add_argument("--kmin",  type=float, default=0.70)
    p.add_argument("--kmax",  type=float, default=1.70)
    p.add_argument("--kstep", type=float, default=0.02)
    args = p.parse_args()

    tdir = Path(args.telescope_dir)
    manifest_path  = tdir / "telescope_manifest.csv"
    patristic_path = tdir / "telescope_patristic.npy"

    assert manifest_path.exists(),  f"Missing {manifest_path}"
    assert patristic_path.exists(), f"Missing {patristic_path}"

    kappa_grid = np.arange(args.kmin, args.kmax + args.kstep / 2, args.kstep)
    kappa_grid = np.round(kappa_grid, 4)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    D_pat_full = np.load(patristic_path)
    print(f"Patristic matrix: {D_pat_full.shape}")

    output = {}

    for label, ckpt_path in [("compact2_k125", args.checkpoint_a),
                              ("redux_k139",   args.checkpoint_b)]:
        print(f"\n{'='*65}")
        print(f"  MODEL: {label}  ({ckpt_path})")
        print(f"{'='*65}")

        model, live_kappa = load_v15model(ckpt_path, device)
        print(f"  live_kappa = {live_kappa:.6f}")

        print(f"\n  Extracting embeddings ...")
        embeddings, valid_idx, trees = extract_embeddings(
            model, device, manifest_path, max_len=args.max_len
        )
        print(f"  {len(valid_idx)} embeddings extracted  (tree dist: {dict(zip(*np.unique(trees, return_counts=True)))})")

        D_pat = D_pat_full[np.ix_(valid_idx, valid_idx)]

        print(f"\n  Dense κ scan ({args.kmin:.2f} → {args.kmax:.2f}, step={args.kstep}) ...")
        scan = dense_scan(embeddings, D_pat, trees, kappa_grid)
        summary = peak_summary(scan, label)

        output[label] = {
            "live_kappa": live_kappa,
            "checkpoint": ckpt_path,
            "n_genomes":  len(valid_idx),
            "summary":    summary,
            "scan":       scan,
        }

        del model
        torch.cuda.empty_cache() if torch.cuda.is_available() else None

    # Side-by-side comparison
    print(f"\n{'='*65}")
    print("  SELF-CONSISTENCY VERDICT")
    print(f"{'='*65}")
    for lbl, res in output.items():
        s = res["summary"]
        print(f"\n  {lbl}  (κ_train={res['live_kappa']:.4f})")
        print(f"    Pearson-log(within) peak: r={s['pearson_log_within_peak']['r']:.4f}"
              f" at κ={s['pearson_log_within_peak']['kappa']:.4f}")
        print(f"    Spearman(within)    peak: ρ={s['spearman_within_peak']['rho']:.4f}"
              f" at κ={s['spearman_within_peak']['kappa']:.4f}")

    print("\n  Hypothesis A (data-intrinsic): peaks agree within ±0.04")
    a_peak = output["compact2_k125"]["summary"]["pearson_log_within_peak"]["kappa"]
    b_peak = output["redux_k139"]["summary"]["pearson_log_within_peak"]["kappa"]
    delta  = abs(a_peak - b_peak)
    print(f"  Peak A = {a_peak:.4f}, Peak B = {b_peak:.4f}, |Δ| = {delta:.4f}")
    print(f"  Verdict: {'Hypothesis A SUPPORTED' if delta < 0.04 else 'Hypothesis B (training bias) possible' if delta > 0.06 else 'Inconclusive'}")

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {args.output}")


if __name__ == "__main__":
    main()
