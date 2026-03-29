#!/usr/bin/env python3
"""
Phase 2 of Telescope Experiment A.

Run ON INFERENCE SERVER. Loads compact2 (V15Model), extracts hyperbolic
embeddings for genomes in the telescope manifest, then fits κ by gradient
descent to minimize stress between Poincaré distances and GTDB patristic
distances.

This is the honest measurement:
  - Reference distances come from GTDB phylogeny (independent of training)
  - Encoder was trained with κ hardcoded to 1.25 (compact2)
  - We ask: does the embedding geometry, evaluated at different curvatures,
    most faithfully represent patristic distances at κ = 1.25?

Usage:
  python3 fit_kappa_telescope.py \\
      --checkpoint /fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt \\
      --telescope-dir ~/telescope/ \\
      --output telescope_kappa_result.json
"""

import argparse, json, math, os, sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from scipy.stats import spearmanr


# ── Poincaré distance (gradient flows through c) ─────────────────────────

def poincare_dist_mat(Z, c, eps=1e-7):
    """Compute full pairwise Poincaré distance matrix."""
    B = Z.shape[0]
    u = Z.unsqueeze(1).expand(B, B, -1)
    v = Z.unsqueeze(0).expand(B, B, -1)
    diff_sq = ((u - v) ** 2).sum(-1)
    u_sq = (u ** 2).sum(-1)
    v_sq = (v ** 2).sum(-1)
    denom = ((1 - c * u_sq) * (1 - c * v_sq)).clamp(min=eps)
    arg = (1 + 2 * c * diff_sq / denom).clamp(min=1 + eps)
    return torch.acosh(arg) / torch.sqrt(c + eps)


# ── Load V15Model ─────────────────────────────────────────────────────────

def load_v15model(checkpoint_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = ckpt.get("model_state_dict", ckpt)

    # Auto-detect architecture from checkpoint
    latent_dim = state["ode_flow.field.0.weight_orig"].shape[1] - 1  # emb_dim - 1
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

    print(f"Checkpoint: latent={latent_dim}, ode_hidden={ode_hidden}, vocab={vocab_size}")
    print(f"Counts: {counts}")

    sys.path.insert(0, str(Path.home()))
    sys.path.insert(0, "/fast/sentrybio/scripts")
    from model_v15_5 import V15Model

    model = V15Model(
        vocab_size=vocab_size,
        latent_dim=latent_dim,
        counts=counts,
        ode_hidden=ode_hidden,
    ).to(device)

    # Remove curvature_history (size mismatch — buffer size depends on training steps)
    state_filtered = {k: v for k, v in state.items() if "curvature_history" not in k}
    missing, unexpected = model.load_state_dict(state_filtered, strict=False)
    if missing:
        print(f"  Missing keys: {len(missing)}")
    if unexpected:
        print(f"  Unexpected keys: {len(unexpected)}")

    model.eval()
    print(f"  live_kappa = {model.live_kappa:.6f}")
    return model, device


# ── Extract embeddings ────────────────────────────────────────────────────

def extract_embeddings(model, device, manifest_path, max_len=8192, seed=42):
    import csv, random
    random.seed(seed)

    rows = []
    with open(manifest_path, newline="") as f:
        for row in csv.DictReader(f):
            tp = row.get("tokenized_path", "").strip()
            if tp and os.path.exists(tp):
                rows.append(row)
            else:
                rows.append(None)  # placeholder: embedding will be NaN

    print(f"Manifest: {len(rows)} entries, "
          f"{sum(r is not None for r in rows)} with existing tokenized files")

    embeddings = []
    valid_idx = []

    with torch.no_grad():
        for i, row in enumerate(rows):
            if row is None:
                continue
            try:
                tokens = np.load(row["tokenized_path"]).astype(np.int64)
                if len(tokens) > max_len:
                    s = random.randint(0, len(tokens) - max_len)
                    tokens = tokens[s:s + max_len]
                t = torch.from_numpy(tokens).unsqueeze(0).to(device)
                # Use encode_angular_only for c-independent geometry
                z_ang = model.encode_angular_only(t)
                embeddings.append(z_ang.squeeze(0).cpu())
                valid_idx.append(i)
            except Exception as e:
                print(f"  Skip {row.get('accession','?')}: {e}")
                continue

    print(f"Extracted {len(embeddings)} embeddings")
    return torch.stack(embeddings), valid_idx


# ── Fit curvature ─────────────────────────────────────────────────────────

def fit_kappa(embeddings, D_pat, c_init_values=None, steps=800, lr=5e-3):
    """
    Find c* = argmin_c Σ_{i<j} (poincare_dist(z_i,z_j;c) - scale * D_pat_{ij})²

    Uses z_ang embeddings (c-independent directions) so that the c-dependence
    comes purely from the distance formula, not the embedding geometry.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Z = embeddings.to(device)
    D = torch.tensor(D_pat, dtype=torch.float32, device=device)

    n = Z.shape[0]
    mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=device), diagonal=1)
    d_ref = D[mask]

    if c_init_values is None:
        c_init_values = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]

    results = []
    for c_init in c_init_values:
        log_c     = torch.tensor(math.log(c_init), requires_grad=True, device=device)
        log_scale = torch.tensor(0.0, requires_grad=True, device=device)

        opt = torch.optim.Adam([log_c, log_scale], lr=lr)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps)

        best_loss = float("inf")
        best_c    = c_init
        for step in range(steps):
            opt.zero_grad()
            c     = torch.exp(log_c)
            scale = torch.exp(log_scale)

            d_hyp  = poincare_dist_mat(Z, c)
            d_flat = d_hyp[mask]
            loss   = F.mse_loss(d_flat, scale * d_ref)
            loss.backward()
            opt.step()
            scheduler.step()

            with torch.no_grad():
                log_c.clamp_(math.log(0.05), math.log(10.0))

            l = loss.item()
            if l < best_loss:
                best_loss = l
                best_c    = torch.exp(log_c).item()

        c_fit     = torch.exp(log_c).item()
        scale_fit = torch.exp(log_scale).item()

        # Spearman correlation
        with torch.no_grad():
            d_fit = poincare_dist_mat(Z, torch.tensor(c_fit, device=device))[mask].cpu().numpy()
        rho, pval = spearmanr(d_fit, d_ref.cpu().numpy())

        results.append({
            "c_init": c_init, "c_fit": c_fit, "c_best": best_c,
            "scale": scale_fit, "stress": best_loss,
            "spearman_rho": float(rho), "p": float(pval),
        })
        print(f"  c_init={c_init:.3f} → c_fit={c_fit:.6f} (best={best_c:.6f})  "
              f"stress={best_loss:.5f}  ρ={rho:.4f}  p={pval:.2e}")

    return results


# ── Evaluate at fixed c values (diagnostic) ──────────────────────────────

def eval_fixed_c(embeddings, D_pat, c_values):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Z = embeddings.to(device)
    D = torch.tensor(D_pat, dtype=torch.float32, device=device)
    n = Z.shape[0]
    mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=device), diagonal=1)
    d_ref_np = D[mask].cpu().numpy()

    rows = []
    for c in c_values:
        with torch.no_grad():
            d_flat = poincare_dist_mat(Z, torch.tensor(c, dtype=torch.float32, device=device))[mask].cpu().numpy()
        rho, pval = spearmanr(d_flat, d_ref_np)
        # Scale-invariant stress: Pearson on log(d+eps)
        from scipy.stats import pearsonr
        pr, pp = pearsonr(np.log(d_flat + 1e-8), np.log(d_ref_np + 1e-8))
        rows.append({"c": c, "spearman_rho": float(rho), "pearson_log": float(pr), "p": float(pval)})
        print(f"  c={c:.4f}  ρ={rho:.4f}  pearson_log={pr:.4f}  p={pval:.2e}")
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",
                   default="/fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt")
    p.add_argument("--telescope-dir", default=os.path.expanduser("~/telescope"))
    p.add_argument("--output",        default="telescope_kappa_result.json")
    p.add_argument("--max-len",  type=int, default=8192)
    p.add_argument("--fit-steps", type=int, default=800)
    p.add_argument("--fit-lr",   type=float, default=5e-3)
    args = p.parse_args()

    tdir = Path(args.telescope_dir)
    manifest_path = tdir / "telescope_manifest.csv"
    patristic_path = tdir / "telescope_patristic.npy"

    assert manifest_path.exists(), f"Missing {manifest_path}"
    assert patristic_path.exists(), f"Missing {patristic_path}"

    print("=" * 65)
    print("TELESCOPE EXPERIMENT A: κ from compact2 + GTDB patristic")
    print("=" * 65)

    # ── Load model ────────────────────────────────────────────────────────
    print("\n[1] Loading compact2 checkpoint...")
    model, device = load_v15model(args.checkpoint)
    training_kappa = model.live_kappa

    # ── Load patristic matrix ─────────────────────────────────────────────
    print(f"\n[2] Loading patristic matrix from {patristic_path}...")
    D_pat_full = np.load(patristic_path)
    print(f"    Shape: {D_pat_full.shape}, range: [{D_pat_full[D_pat_full>0].min():.4f}, {D_pat_full.max():.4f}]")

    # ── Extract embeddings ────────────────────────────────────────────────
    print(f"\n[3] Extracting embeddings (max_len={args.max_len})...")
    embeddings, valid_idx = extract_embeddings(model, device, manifest_path, max_len=args.max_len)

    # Subset patristic matrix to valid indices only
    D_pat = D_pat_full[np.ix_(valid_idx, valid_idx)]
    n = len(valid_idx)
    print(f"    Using {n}×{n} submatrix of patristic distances")

    # ── Diagnostic: evaluate at fixed κ values ───────────────────────────
    print(f"\n[4] Evaluating at fixed κ values (embedding quality probe)...")
    scan_c = [0.5, 0.75, 0.88, 1.0, 1.1, 1.25, 1.5, 2.0]
    fixed_evals = eval_fixed_c(embeddings, D_pat, scan_c)

    # ── Fit κ ─────────────────────────────────────────────────────────────
    print(f"\n[5] Fitting κ from {len([0.5,0.75,1.0,1.25,1.5,2.0])} initializations "
          f"({args.fit_steps} steps each)...")
    fit_results = fit_kappa(embeddings, D_pat, steps=args.fit_steps, lr=args.fit_lr)

    # ── Summary ───────────────────────────────────────────────────────────
    c_fits = [r["c_fit"] for r in fit_results]
    c_bests = [r["c_best"] for r in fit_results]
    theory = (1.6 * math.log(2)) ** 2

    rhos = [r["spearman_rho"] for r in fit_results]
    best_run = max(fit_results, key=lambda r: r["spearman_rho"])

    print("\n" + "=" * 65)
    print("RESULTS")
    print("=" * 65)
    print(f"N genomes:              {n}")
    print(f"κ (training, compact2): {training_kappa:.6f}")
    print(f"κ theory (h·ln2)²:      {theory:.6f}")
    print(f"κ fit mean:             {np.mean(c_fits):.6f} ± {np.std(c_fits):.6f}")
    print(f"κ fit median:           {np.median(c_fits):.6f}")
    print(f"κ best (lowest stress): {min(c_bests, key=lambda c: abs(c - np.mean(c_bests))):.6f}")
    print(f"Best Spearman ρ:        {best_run['spearman_rho']:.4f}  (c={best_run['c_fit']:.4f})")
    print(f"Agreement with theory:  {abs(np.mean(c_fits)-theory)/theory*100:.1f}%")

    # Fixed-c Spearman peak
    best_fixed = max(fixed_evals, key=lambda r: r["spearman_rho"])
    print(f"\nFixed-κ scan peak:      ρ={best_fixed['spearman_rho']:.4f} at κ={best_fixed['c']:.4f}")

    out = {
        "n_genomes": n,
        "training_kappa": training_kappa,
        "theory_kappa": theory,
        "kappa_fit_mean": float(np.mean(c_fits)),
        "kappa_fit_std":  float(np.std(c_fits)),
        "kappa_fit_median": float(np.median(c_fits)),
        "best_spearman_rho": float(best_run["spearman_rho"]),
        "best_spearman_c": float(best_run["c_fit"]),
        "fixed_kappa_scan": fixed_evals,
        "fit_runs": fit_results,
    }
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {args.output}")


if __name__ == "__main__":
    main()
