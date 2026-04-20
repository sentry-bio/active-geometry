#!/usr/bin/env python3
"""
Eukaryote Telescope Experiment.

Sweeps κ over [0.7, 2.0] and measures Spearman ρ between:
  - Poincaré distances (compact2 embeddings at each κ)
  - OG profile distances (eggNOG OG cosine distance as molecular proxy)

This is the eukaryote-only version of Telescope Experiment A (which used
GTDB patristic distances for prokaryotes). Since no published eukaryotic
ML tree with branch lengths is available, we use eggNOG scale3 OG profile
cosine distance as a phylogenetically-informed molecular distance proxy.

Reference: eggNOG family-level (scale3) OG profiles, 200 OGs.
Encoder: compact2 (V15Model, best.pt, trained at κ=1.2505).

Usage:
  python3 euk_telescope.py \\
      --checkpoint /fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt \\
      --og-profiles /home/rohit/eggnog/og_level_profiles.npz \\
      --manifest /home/rohit/eggnog/manifest_go_kegg_v2.csv \\
      --output euk_telescope_result.json \\
      --n-genomes 300
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from scipy.stats import spearmanr, pearsonr


# ── Poincaré distance ────────────────────────────────────────────────────────

def poincare_dist_mat(Z, c, eps=1e-7):
    """Full pairwise Poincaré distance matrix. c = κ."""
    B = Z.shape[0]
    u = Z.unsqueeze(1).expand(B, B, -1)
    v = Z.unsqueeze(0).expand(B, B, -1)
    diff_sq = ((u - v) ** 2).sum(-1)
    u_sq = (u ** 2).sum(-1)
    v_sq = (v ** 2).sum(-1)
    denom = ((1 - c * u_sq) * (1 - c * v_sq)).clamp(min=eps)
    arg = (1 + 2 * c * diff_sq / denom).clamp(min=1 + eps)
    return torch.acosh(arg) / torch.sqrt(c + eps)


# ── Load V15Model ─────────────────────────────────────────────────────────────

def load_v15model(checkpoint_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # GPU safety: fall back to CPU if less than 25% free
    if device.type == "cuda":
        total = torch.cuda.get_device_properties(0).total_memory
        reserved = torch.cuda.memory_reserved(0)
        free_frac = (total - reserved) / total
        if free_frac < 0.25:
            print(f"  GPU only {free_frac*100:.0f}% free — using CPU")
            device = torch.device("cpu")

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

    print(f"  latent={latent_dim}, ode_hidden={ode_hidden}, vocab={vocab_size}")
    print(f"  prototypes: {counts}")

    sys.path.insert(0, "/home/rohit")
    sys.path.insert(0, "/fast/sentrybio/biosphere-atlas/biosphere_atlas/core")
    from model_v15_5 import V15Model

    model = V15Model(
        vocab_size=vocab_size,
        latent_dim=latent_dim,
        counts=counts,
        ode_hidden=ode_hidden,
    ).to(device)

    state_filtered = {k: v for k, v in state.items() if "curvature_history" not in k}
    missing, unexpected = model.load_state_dict(state_filtered, strict=False)
    if missing:
        print(f"  Missing keys: {len(missing)}")
    if unexpected:
        print(f"  Unexpected keys: {len(unexpected)}")

    model.eval()
    print(f"  live_kappa = {model.live_kappa:.6f}")
    return model, device


# ── Select diverse eukaryote sample ──────────────────────────────────────────

def select_diverse_sample(manifest_path, og_accession_set, n_target, seed=42):
    """
    Stratified sample across phyla, all with tokenized files and OG profiles.
    Returns list of dicts with accession, tokenized_path, phylum, species, og_idx.
    """
    import pandas as pd
    import random
    random.seed(seed)
    np.random.seed(seed)

    df = pd.read_csv(manifest_path,
                     usecols=['accession', 'domain', 'tokenized_path',
                               'species', 'phylum'],
                     low_memory=False)

    euk = df[df['domain'] == 'Eukaryota'].copy()
    euk = euk[euk['accession'].isin(og_accession_set)]
    euk = euk[euk['tokenized_path'].apply(lambda p: bool(p) and os.path.exists(str(p)))]
    euk = euk.reset_index(drop=True)

    print(f"  Eligible eukaryotes: {len(euk)}")
    print(f"  Phylum distribution:")
    print(euk['phylum'].value_counts().head(10).to_string())

    # Stratified sample: proportional to phylum, min 5 per phylum
    phyla = euk['phylum'].value_counts()
    sample_rows = []
    for phylum, count in phyla.items():
        phylum_df = euk[euk['phylum'] == phylum]
        n_phylum = max(5, int(round(n_target * count / len(euk))))
        n_phylum = min(n_phylum, len(phylum_df))
        sample_rows.append(phylum_df.sample(n_phylum, random_state=seed))

    sample = pd.concat(sample_rows).drop_duplicates(subset='accession')
    # Trim if over target
    if len(sample) > n_target:
        sample = sample.sample(n_target, random_state=seed)

    sample = sample.reset_index(drop=True)
    print(f"\n  Selected {len(sample)} genomes across {sample['phylum'].nunique()} phyla")
    print(sample['phylum'].value_counts().head(10).to_string())

    return sample


# ── Extract embeddings ────────────────────────────────────────────────────────

def extract_embeddings(model, device, sample_df, max_len=8192, seed=42):
    import random
    random.seed(seed)

    embeddings = []
    valid_mask = []

    with torch.no_grad():
        for i, row in sample_df.iterrows():
            try:
                tokens = np.load(row['tokenized_path']).astype(np.int64)
                if len(tokens) > max_len:
                    s = random.randint(0, len(tokens) - max_len)
                    tokens = tokens[s:s + max_len]
                t = torch.from_numpy(tokens).unsqueeze(0).to(device)
                z = model.encode_angular_only(t)
                embeddings.append(z.squeeze(0).cpu())
                valid_mask.append(True)
            except Exception as e:
                print(f"  Skip {row.get('accession','?')}: {e}")
                valid_mask.append(False)

    n_ok = sum(valid_mask)
    print(f"  Extracted {n_ok}/{len(sample_df)} embeddings")
    return torch.stack(embeddings), np.array(valid_mask)


# ── Build OG profile distance matrix ─────────────────────────────────────────

def build_og_distance_matrix(og_data, sample_accessions, scale='scale3'):
    """
    OG profile cosine distance: 1 - cosine_similarity.
    scale3 = family-level OGs, best phylogenetic resolution for eukaryotes.
    """
    all_acc = list(og_data['accessions'])
    acc_idx = {a: i for i, a in enumerate(all_acc)}

    profiles_key = f"{scale}_profiles"
    profiles = og_data[profiles_key]  # (N_all, 200)

    # Get row indices for our sample
    idx = []
    missing = []
    for acc in sample_accessions:
        if acc in acc_idx:
            idx.append(acc_idx[acc])
        else:
            idx.append(None)
            missing.append(acc)

    if missing:
        print(f"  WARNING: {len(missing)} accessions not in OG profiles")

    # Replace missing with mean profile
    mean_prof = profiles.mean(axis=0)
    P = np.stack([profiles[i] if i is not None else mean_prof for i in idx])

    # L2-normalize for cosine distance
    norms = np.linalg.norm(P, axis=1, keepdims=True).clip(min=1e-8)
    P_norm = P / norms

    # Cosine distance = 1 - cosine_similarity
    cos_sim = P_norm @ P_norm.T
    D = (1.0 - cos_sim).clip(min=0.0)

    # Force diagonal to zero
    np.fill_diagonal(D, 0.0)

    dist_range = D[D > 0]
    print(f"  OG {scale} distance: min={dist_range.min():.4f}, "
          f"max={D.max():.4f}, mean={dist_range.mean():.4f}")

    return D


# ── Evaluate at fixed κ values ────────────────────────────────────────────────

def eval_fixed_kappa(embeddings, D_og, kappa_values):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Z = embeddings.to(device)
    D = torch.tensor(D_og, dtype=torch.float32, device=device)
    n = Z.shape[0]
    mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=device), diagonal=1)
    d_ref_np = D[mask].cpu().numpy()

    rows = []
    for kappa in kappa_values:
        c = torch.tensor(kappa, dtype=torch.float32, device=device)
        with torch.no_grad():
            d_hyp = poincare_dist_mat(Z, c)[mask].cpu().numpy()
        rho, pval = spearmanr(d_hyp, d_ref_np)
        pr, _ = pearsonr(np.log(d_hyp + 1e-8), np.log(d_ref_np + 1e-8))
        rows.append({
            "kappa": float(kappa),
            "spearman_rho": float(rho),
            "pearson_log": float(pr),
            "p": float(pval),
        })
        print(f"  κ={kappa:.3f}  ρ={rho:.4f}  pearson_log={pr:.4f}  p={pval:.2e}")
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",
                   default="/fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt")
    p.add_argument("--og-profiles",
                   default="/home/rohit/eggnog/og_level_profiles.npz")
    p.add_argument("--manifest",
                   default="/home/rohit/eggnog/manifest_go_kegg_v2.csv")
    p.add_argument("--output", default="euk_telescope_result.json")
    p.add_argument("--n-genomes", type=int, default=300)
    p.add_argument("--max-len",   type=int, default=8192)
    p.add_argument("--og-scale",  default="scale3",
                   choices=["scale1", "scale2", "scale3"])
    p.add_argument("--seed",      type=int, default=42)
    args = p.parse_args()

    print("=" * 65)
    print("EUKARYOTE TELESCOPE: κ from compact2 + eggNOG OG distances")
    print("=" * 65)

    # ── Load OG profiles ───────────────────────────────────────────────────
    print(f"\n[1] Loading OG profiles from {args.og_profiles}...")
    og_data = np.load(args.og_profiles, allow_pickle=True)
    og_accessions = set(og_data['accessions'])
    print(f"    {len(og_accessions)} accessions in OG profile database")

    # ── Select sample ──────────────────────────────────────────────────────
    print(f"\n[2] Selecting {args.n_genomes} diverse eukaryotes...")
    sample_df = select_diverse_sample(
        args.manifest, og_accessions, args.n_genomes, args.seed
    )

    # ── Load model ─────────────────────────────────────────────────────────
    print(f"\n[3] Loading compact2 checkpoint...")
    model, device = load_v15model(args.checkpoint)
    training_kappa = model.live_kappa

    # ── Extract embeddings ─────────────────────────────────────────────────
    print(f"\n[4] Extracting embeddings (max_len={args.max_len})...")
    embeddings, valid_mask = extract_embeddings(
        model, device, sample_df, max_len=args.max_len, seed=args.seed
    )
    valid_df = sample_df[valid_mask].reset_index(drop=True)
    valid_accessions = valid_df['accession'].tolist()
    n = len(valid_accessions)
    print(f"    {n} valid embeddings")

    # ── Build OG distance matrix ───────────────────────────────────────────
    print(f"\n[5] Building OG profile distance matrix (scale={args.og_scale})...")
    D_og = build_og_distance_matrix(og_data, valid_accessions, scale=args.og_scale)

    # Also build scale1 and scale2 for comparison
    print("    Also building scale1 (universal) and scale2 (domain) for comparison...")
    D_og_s1 = build_og_distance_matrix(og_data, valid_accessions, scale="scale1")
    D_og_s2 = build_og_distance_matrix(og_data, valid_accessions, scale="scale2")

    # ── Coarse κ sweep ─────────────────────────────────────────────────────
    print(f"\n[6] Coarse κ sweep [0.70, 2.00] using {args.og_scale} OG distances...")
    coarse_kappas = [0.70, 0.80, 0.88, 1.00, 1.10, 1.25, 1.35, 1.40, 1.50, 2.00]
    coarse_results = eval_fixed_kappa(embeddings, D_og, coarse_kappas)

    # Find coarse peak
    peak_coarse = max(coarse_results, key=lambda r: r["spearman_rho"])
    peak_kappa = peak_coarse["kappa"]
    print(f"\n  Coarse peak: κ={peak_kappa:.2f}, ρ={peak_coarse['spearman_rho']:.4f}")

    # ── Fine κ sweep around peak ────────────────────────────────────────────
    fine_lo = max(0.7, peak_kappa - 0.15)
    fine_hi = min(2.0, peak_kappa + 0.15)
    fine_kappas = list(np.arange(fine_lo, fine_hi + 0.01, 0.02))
    print(f"\n[7] Fine κ sweep [{fine_lo:.2f}, {fine_hi:.2f}] using {args.og_scale} OG distances...")
    fine_results = eval_fixed_kappa(embeddings, D_og, fine_kappas)

    peak_fine = max(fine_results, key=lambda r: r["spearman_rho"])
    peak_kappa_fine = peak_fine["kappa"]
    print(f"\n  Fine peak: κ={peak_kappa_fine:.3f}, ρ={peak_fine['spearman_rho']:.4f}")

    # ── Compare across OG scales ───────────────────────────────────────────
    print(f"\n[8] Comparing κ sweep across all 3 OG scales (at coarse grid)...")
    print("  scale1 (universal):")
    s1_results = eval_fixed_kappa(embeddings, D_og_s1, coarse_kappas)
    print("  scale2 (domain):")
    s2_results = eval_fixed_kappa(embeddings, D_og_s2, coarse_kappas)

    # ── Summary ────────────────────────────────────────────────────────────
    theory = (1.6 * math.log(2)) ** 2
    theory_h170 = (1.70 * math.log(2)) ** 2

    print("\n" + "=" * 65)
    print("RESULTS")
    print("=" * 65)
    print(f"N eukaryote genomes:      {n}")
    print(f"OG reference scale:       {args.og_scale}")
    print(f"κ (compact2 training):    {training_kappa:.6f}")
    print(f"κ theory (h=1.60):        {theory:.6f}")
    print(f"κ theory (h=1.70):        {theory_h170:.6f}")
    print(f"κ peak (coarse):          {peak_kappa:.3f}  (ρ={peak_coarse['spearman_rho']:.4f})")
    print(f"κ peak (fine):            {peak_kappa_fine:.3f}  (ρ={peak_fine['spearman_rho']:.4f})")

    # Phylum breakdown of sample
    phylum_counts = valid_df['phylum'].value_counts().to_dict()

    out = {
        "experiment": "euk_telescope",
        "date": "2026-03-12",
        "n_genomes": n,
        "n_requested": args.n_genomes,
        "og_scale": args.og_scale,
        "phylum_composition": phylum_counts,
        "training_kappa": float(training_kappa),
        "theory_kappa_h160": float(theory),
        "theory_kappa_h170": float(theory_h170),
        "coarse_scan": {
            "kappa":        [r["kappa"] for r in coarse_results],
            "spearman":     [r["spearman_rho"] for r in coarse_results],
            "pearson_log":  [r["pearson_log"] for r in coarse_results],
        },
        "fine_scan": {
            "kappa":        [r["kappa"] for r in fine_results],
            "spearman":     [r["spearman_rho"] for r in fine_results],
            "pearson_log":  [r["pearson_log"] for r in fine_results],
        },
        "peak_kappa_coarse": float(peak_kappa),
        "peak_rho_coarse": float(peak_coarse["spearman_rho"]),
        "peak_kappa_fine": float(peak_kappa_fine),
        "peak_rho_fine": float(peak_fine["spearman_rho"]),
        "scale_comparison": {
            "scale1": {
                "kappa":    [r["kappa"] for r in s1_results],
                "spearman": [r["spearman_rho"] for r in s1_results],
            },
            "scale2": {
                "kappa":    [r["kappa"] for r in s2_results],
                "spearman": [r["spearman_rho"] for r in s2_results],
            },
            "scale3": {
                "kappa":    [r["kappa"] for r in coarse_results],
                "spearman": [r["spearman_rho"] for r in coarse_results],
            },
        },
    }

    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {args.output}")


if __name__ == "__main__":
    main()
