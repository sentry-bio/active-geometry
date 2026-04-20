#!/usr/bin/env python3
"""
Fungal Telescope Experiment — canonical eukaryote κ measurement.

Sweeps κ over [0.7, 2.0] and measures Spearman ρ between:
  - Poincaré distances from compact2 embeddings at each κ
  - Patristic distances from Li et al. 2021 fungal ML tree
    (1,672 taxa, 290 BUSCO genes, IQ-TREE ultrafast bootstrap)

This is the eukaryote analog of Telescope Experiment A (GTDB prokaryotes).
Reference tree has real ML-estimated branch lengths — not a proxy.

Usage:
  python3 fungal_telescope.py \\
      --checkpoint /fast/sentrybio/checkpoints/v15_5_compact2_diagnostic/best.pt \\
      --tree /home/rohit/li2021_fungal_tree.treefile \\
      --manifest /home/rohit/eggnog/manifest_go_kegg_v2.csv \\
      --output /home/rohit/fungal_telescope_result.json
"""

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from scipy.stats import spearmanr, pearsonr


# ── Poincaré distance ─────────────────────────────────────────────────────────

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


# ── Load V15Model ─────────────────────────────────────────────────────────────

def load_v15model(checkpoint_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type == "cuda":
        total = torch.cuda.get_device_properties(0).total_memory
        reserved = torch.cuda.memory_reserved(0)
        free_frac = (total - reserved) / total
        if free_frac < 0.20:
            print(f"  GPU only {free_frac*100:.0f}% free — using CPU")
            device = torch.device("cpu")

    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = ckpt.get("model_state_dict", ckpt)

    latent_dim = state["ode_flow.field.0.weight_orig"].shape[1] - 1
    ode_hidden  = state["ode_flow.field.0.bias"].shape[0]
    counts = {k: state[f"{k}.prototypes"].shape[0]
              for k in ["bact_fam","arch_fam","euk_fam","bact_gen","arch_gen","euk_gen"]}
    vocab_size = state["encoder.encoder.embed.weight"].shape[0]

    print(f"  latent={latent_dim}, ode_hidden={ode_hidden}, vocab={vocab_size}")

    sys.path.insert(0, "/home/rohit/biosphere_inference")
    sys.path.insert(0, "/home/rohit")
    from model_v15_5 import V15Model

    model = V15Model(
        vocab_size=vocab_size, latent_dim=latent_dim,
        counts=counts, ode_hidden=ode_hidden,
    ).to(device)
    state_filtered = {k: v for k, v in state.items() if "curvature_history" not in k}
    missing, unexpected = model.load_state_dict(state_filtered, strict=False)
    if missing:    print(f"  Missing keys: {len(missing)}")
    if unexpected: print(f"  Unexpected keys: {len(unexpected)}")
    model.eval()
    print(f"  live_kappa = {model.live_kappa:.6f}")
    return model, device


# ── Build patristic distance matrix from ML tree ──────────────────────────────

def build_patristic_matrix(tree_path, target_species_list):
    """
    Compute pairwise patristic distances (sum of branch lengths on path between
    each pair) for the species in target_species_list.

    Uses DendroPy for tree I/O and PDM computation.
    Returns: (D, valid_indices) where D is [n,n] numpy array, valid_indices
    maps row/col back to positions in target_species_list.
    """
    import dendropy

    print(f"  Loading tree from {tree_path}...")
    tree = dendropy.Tree.get(path=tree_path, schema="newick",
                              preserve_underscores=True)

    # Build label → taxon map
    taxon_labels = {t.label: t for t in tree.taxon_namespace}
    print(f"  Tree taxa: {len(taxon_labels)}")

    # Normalize tree labels: underscore→space, lowercase
    taxon_norm = {lbl.replace('_', ' ').lower(): lbl for lbl in taxon_labels}

    # Match target species
    matched = []   # (idx_in_target, tree_label)
    unmatched = []
    for i, sp in enumerate(target_species_list):
        sp_norm = sp.lower().strip()
        if sp_norm in taxon_norm:
            matched.append((i, taxon_norm[sp_norm]))
        else:
            unmatched.append(sp)

    print(f"  Matched: {len(matched)}/{len(target_species_list)}")
    if unmatched[:5]:
        print(f"  Unmatched examples: {unmatched[:5]}")

    if len(matched) < 10:
        raise ValueError("Too few species matched in tree (<10). Check species names.")

    # Compute pairwise patristic distances using DendroPy PDM
    print(f"  Computing {len(matched)}×{len(matched)} patristic distance matrix...")
    pdm = tree.phylogenetic_distance_matrix()

    n = len(matched)
    D = np.zeros((n, n), dtype=np.float32)
    for i, (_, lbl_i) in enumerate(matched):
        tax_i = taxon_labels[lbl_i]
        for j, (_, lbl_j) in enumerate(matched):
            if i == j:
                continue
            tax_j = taxon_labels[lbl_j]
            D[i, j] = pdm.distance(tax_i, tax_j)

    valid_target_indices = [idx for idx, _ in matched]
    matched_labels = [lbl for _, lbl in matched]

    d_vals = D[D > 0]
    print(f"  Patristic distance range: [{d_vals.min():.4f}, {D.max():.4f}], "
          f"mean={d_vals.mean():.4f}")
    return D, valid_target_indices, matched_labels


# ── Extract embeddings ────────────────────────────────────────────────────────

def extract_embeddings(model, device, tok_paths, max_len=8192, seed=42):
    import random
    random.seed(seed)

    embeddings = []
    valid_mask = []

    with torch.no_grad():
        for i, tp in enumerate(tok_paths):
            try:
                tokens = np.load(tp).astype(np.int64)
                if len(tokens) > max_len:
                    s = random.randint(0, len(tokens) - max_len)
                    tokens = tokens[s:s + max_len]
                t = torch.from_numpy(tokens).unsqueeze(0).to(device)
                z = model.encode_angular_only(t)
                embeddings.append(z.squeeze(0).cpu())
                valid_mask.append(True)
            except Exception as e:
                print(f"  Skip [{i}]: {e}")
                valid_mask.append(False)

            if (i + 1) % 100 == 0:
                print(f"  [{i+1}/{len(tok_paths)}] embedded")

    n_ok = sum(valid_mask)
    print(f"  Extracted {n_ok}/{len(tok_paths)} embeddings")
    return torch.stack(embeddings), np.array(valid_mask)


# ── Evaluate at fixed κ values ────────────────────────────────────────────────

def eval_fixed_kappa(embeddings, D_pat, kappa_values):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Z = embeddings.to(device)
    D = torch.tensor(D_pat, dtype=torch.float32, device=device)
    n = Z.shape[0]
    mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=device), diagonal=1)
    d_ref_np = D[mask].cpu().numpy()

    rows = []
    for kappa in kappa_values:
        c = torch.tensor(kappa, dtype=torch.float32, device=device)
        with torch.no_grad():
            d_hyp = poincare_dist_mat(Z, c)[mask].cpu().numpy()
        rho, pval = spearmanr(d_hyp, d_ref_np)
        # Log-space Pearson (scale-invariant)
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
    p.add_argument("--tree",     default="/home/rohit/li2021_fungal_tree.treefile")
    p.add_argument("--manifest", default="/home/rohit/eggnog/manifest_go_kegg_v2.csv")
    p.add_argument("--output",   default="/home/rohit/fungal_telescope_result.json")
    p.add_argument("--max-len",  type=int, default=8192)
    p.add_argument("--max-genomes", type=int, default=0,
                   help="Cap genomes for testing (0 = use all matched)")
    p.add_argument("--seed",     type=int, default=42)
    args = p.parse_args()

    print("=" * 65)
    print("FUNGAL TELESCOPE: κ from compact2 + Li2021 ML patristic distances")
    print("=" * 65)

    import pandas as pd

    # ── Load manifest (fungi only) ─────────────────────────────────────────
    print("\n[1] Loading fungal genomes from manifest...")
    df = pd.read_csv(args.manifest,
                     usecols=['accession','domain','species','phylum','tokenized_path'],
                     low_memory=False)
    fungi = df[
        (df['domain'] == 'Eukaryota') &
        (df['phylum'].isin(['Ascomycota','Basidiomycota','Mucoromycota','Microsporidia']))
    ].copy()
    fungi = fungi[fungi['tokenized_path'].apply(
        lambda p: bool(p) and isinstance(p, str) and os.path.exists(p)
    )].reset_index(drop=True)
    print(f"  Fungal genomes with tokenized files: {len(fungi)}")
    print(f"  Phyla: {fungi['phylum'].value_counts().to_dict()}")

    # ── Build patristic matrix ─────────────────────────────────────────────
    print(f"\n[2] Loading ML tree and building patristic distance matrix...")
    D_pat, valid_target_idx, matched_labels = build_patristic_matrix(
        args.tree, fungi['species'].fillna('').tolist()
    )

    # Subset manifest to matched species
    matched_df = fungi.iloc[valid_target_idx].reset_index(drop=True)
    matched_df['tree_label'] = matched_labels

    if args.max_genomes > 0 and len(matched_df) > args.max_genomes:
        import random
        random.seed(args.seed)
        keep = sorted(random.sample(range(len(matched_df)), args.max_genomes))
        matched_df = matched_df.iloc[keep].reset_index(drop=True)
        D_pat = D_pat[np.ix_(keep, keep)]
        print(f"  Subsampled to {args.max_genomes} genomes")

    n = len(matched_df)
    print(f"  Final corpus: {n} genomes")
    print(f"  Phylum breakdown: {matched_df['phylum'].value_counts().to_dict()}")

    # Normalize patristic distances by max (same as prokaryote telescope)
    d_max = D_pat.max()
    D_norm = D_pat / d_max
    print(f"  Patristic max={d_max:.4f}, after norm max=1.000")

    # ── Load model ─────────────────────────────────────────────────────────
    print(f"\n[3] Loading compact2 checkpoint...")
    model, device = load_v15model(args.checkpoint)
    training_kappa = model.live_kappa

    # ── Extract embeddings ─────────────────────────────────────────────────
    print(f"\n[4] Extracting embeddings ({n} genomes, max_len={args.max_len})...")
    tok_paths = matched_df['tokenized_path'].tolist()
    embeddings, valid_mask = extract_embeddings(
        model, device, tok_paths, max_len=args.max_len, seed=args.seed
    )

    # Filter to successfully embedded
    if not valid_mask.all():
        valid_idx = np.where(valid_mask)[0]
        embeddings = embeddings  # already stacked from valid ones only
        D_norm = D_norm[np.ix_(valid_idx, valid_idx)]
        matched_df = matched_df.iloc[valid_idx].reset_index(drop=True)
        n = len(matched_df)
        print(f"  After filtering: {n} valid genomes")

    # ── Coarse κ sweep ─────────────────────────────────────────────────────
    print(f"\n[5] Coarse κ sweep [0.70, 2.00]...")
    coarse_kappas = [0.70, 0.80, 0.88, 1.00, 1.10, 1.25, 1.35, 1.40, 1.50, 2.00]
    coarse_results = eval_fixed_kappa(embeddings, D_norm, coarse_kappas)

    peak_coarse = max(coarse_results, key=lambda r: r["spearman_rho"])
    peak_kappa = peak_coarse["kappa"]
    print(f"\n  Coarse peak: κ={peak_kappa:.2f}, ρ={peak_coarse['spearman_rho']:.4f}")

    # ── Fine κ sweep ───────────────────────────────────────────────────────
    fine_lo = max(0.70, peak_kappa - 0.15)
    fine_hi = min(2.00, peak_kappa + 0.15)
    fine_kappas = list(np.arange(fine_lo, fine_hi + 0.01, 0.02))
    print(f"\n[6] Fine κ sweep [{fine_lo:.2f}, {fine_hi:.2f}]...")
    fine_results = eval_fixed_kappa(embeddings, D_norm, fine_kappas)

    peak_fine = max(fine_results, key=lambda r: r["spearman_rho"])
    peak_kappa_fine = peak_fine["kappa"]
    print(f"\n  Fine peak: κ={peak_kappa_fine:.3f}, ρ={peak_fine['spearman_rho']:.4f}")

    # ── Summary ────────────────────────────────────────────────────────────
    theory_160 = (1.60 * math.log(2)) ** 2
    theory_170 = (1.70 * math.log(2)) ** 2

    print("\n" + "=" * 65)
    print("RESULTS")
    print("=" * 65)
    print(f"N fungal genomes:          {n}")
    print(f"κ (compact2 training):     {training_kappa:.6f}")
    print(f"κ theory (h=1.60·ln2)²:   {theory_160:.6f}")
    print(f"κ theory (h=1.70·ln2)²:   {theory_170:.6f}")
    print(f"κ GTDB prokaryote peak:    ~1.34–1.39 (Telescope A)")
    print(f"κ peak coarse:             {peak_kappa:.3f}  (ρ={peak_coarse['spearman_rho']:.4f})")
    print(f"κ peak fine:               {peak_kappa_fine:.3f}  (ρ={peak_fine['spearman_rho']:.4f})")

    phylum_counts = matched_df['phylum'].value_counts().to_dict()

    out = {
        "experiment": "fungal_telescope",
        "reference": "Li et al. 2021 Current Biology — 1,672 taxa, 290 BUSCO genes, IQ-TREE",
        "date": "2026-03-13",
        "n_genomes": n,
        "phylum_composition": phylum_counts,
        "training_kappa": float(training_kappa),
        "theory_kappa_h160": float(theory_160),
        "theory_kappa_h170": float(theory_170),
        "patristic_max_raw": float(d_max),
        "coarse_scan": {
            "kappa":       [r["kappa"] for r in coarse_results],
            "spearman":    [r["spearman_rho"] for r in coarse_results],
            "pearson_log": [r["pearson_log"] for r in coarse_results],
        },
        "fine_scan": {
            "kappa":       [r["kappa"] for r in fine_results],
            "spearman":    [r["spearman_rho"] for r in fine_results],
            "pearson_log": [r["pearson_log"] for r in fine_results],
        },
        "peak_kappa_coarse": float(peak_kappa),
        "peak_rho_coarse": float(peak_coarse["spearman_rho"]),
        "peak_kappa_fine": float(peak_kappa_fine),
        "peak_rho_fine": float(peak_fine["spearman_rho"]),
    }

    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {args.output}")


if __name__ == "__main__":
    main()
