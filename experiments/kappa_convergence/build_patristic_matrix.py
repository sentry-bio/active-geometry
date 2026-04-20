#!/usr/bin/env python3
"""
Phase 1 of Telescope Experiment A.

Run LOCALLY. Loads GTDB bac120 + ar53 trees, finds overlap with tokenized
genomes on the inference server, samples N diverse leaves, computes pairwise
patristic distances, and saves:
  - telescope_manifest.csv   : accession, gtdb_leaf, tokenized_path, domain, phylum, genus
  - telescope_patristic.npy  : float32 [N x N] patristic distance matrix
  - telescope_leaves.txt     : ordered list of accessions (row/col index)
"""

import argparse, csv, gzip, os, random, re, sys
from io import StringIO
from pathlib import Path

import numpy as np
from Bio import Phylo


GTDB_DIR = Path("/Users/rohitfenn/wolframaplpha/biosphere-geometry-of-life-main copy/data/trees/gtdb")
TOKENIZED_LIST = "/tmp/tokenized_list.txt"
TOKENIZED_BASE = "/fast/sentrybio/data/tokenized"  # path on inference server


def load_gtdb_leaves(tree_path):
    """Return dict: accession → (gtdb_leaf, branch_lengths_tree_obj)."""
    print(f"Loading {tree_path.name} ...", end=" ", flush=True)
    with gzip.open(tree_path, "rt") as f:
        content = f.read()
    # Extract leaf names via regex (fast, avoids full parse first)
    leaves_raw = re.findall(r"((?:RS|GB)_GC[AF]_\d+\.\d+)", content)
    lookup = {leaf.split("_", 1)[1]: leaf for leaf in set(leaves_raw)}
    print(f"{len(lookup)} leaves")
    return lookup, content


def parse_tree(newick_str):
    return Phylo.read(StringIO(newick_str), "newick")


def compute_patristic_distances(tree, leaves, leaf_to_id):
    """
    Compute all-pairs patristic distances for the given leaf names.
    Returns float32 [N x N] matrix.
    """
    n = len(leaves)
    D = np.zeros((n, n), dtype=np.float32)

    # Build {name: Clade} mapping from the tree
    name_to_clade = {c.name: c for c in tree.find_clades() if c.name}

    print(f"  Computing {n*(n-1)//2} pairwise patristic distances...", flush=True)
    for i in range(n):
        if i % 20 == 0:
            print(f"    row {i}/{n}", flush=True)
        ci = name_to_clade.get(leaves[i])
        if ci is None:
            continue
        for j in range(i + 1, n):
            cj = name_to_clade.get(leaves[j])
            if cj is None:
                continue
            d = tree.distance(ci, cj)
            D[i, j] = D[j, i] = d

    return D


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-samples", type=int, default=300,
                   help="Target number of genomes to sample")
    p.add_argument("--bac-frac", type=float, default=0.7,
                   help="Fraction from bacteria (rest from archaea)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-dir", default="/tmp/telescope")
    p.add_argument("--manifest", default=None,
                   help="Path to manifest CSV on inference (for taxonomy labels). "
                        "Leave None to skip taxonomy columns.")
    args = p.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load tokenized list ────────────────────────────────────────────────
    with open(TOKENIZED_LIST) as f:
        tokenized = set(l.strip() for l in f)
    print(f"Tokenized genomes on inference: {len(tokenized)}")

    # ── Load GTDB trees ────────────────────────────────────────────────────
    bac_lookup, bac_nwk = load_gtdb_leaves(GTDB_DIR / "bac120.tree.gz")
    ar_lookup, ar_nwk  = load_gtdb_leaves(GTDB_DIR / "ar53.tree.gz")

    bac_overlap = sorted(tokenized & set(bac_lookup.keys()))
    ar_overlap  = sorted(tokenized & set(ar_lookup.keys()))
    print(f"Bacterial overlap: {len(bac_overlap)}")
    print(f"Archaeal  overlap: {len(ar_overlap)}")

    # ── Sample ────────────────────────────────────────────────────────────
    n_bac = min(int(args.n_samples * args.bac_frac), len(bac_overlap))
    n_ar  = min(args.n_samples - n_bac, len(ar_overlap))
    sam_bac = random.sample(bac_overlap, n_bac)
    sam_ar  = random.sample(ar_overlap, n_ar)
    print(f"Sampled: {n_bac} bacteria + {n_ar} archaea = {n_bac+n_ar} total")

    # Build ordered list
    # For bacteria: leaf name in bac120 tree
    # For archaea:  leaf name in ar53 tree
    entries = []
    for acc in sam_bac:
        entries.append({"accession": acc, "gtdb_leaf": bac_lookup[acc],
                        "tree": "bac120", "tokenized_path":
                        f"{TOKENIZED_BASE}/{acc}.npy"})
    for acc in sam_ar:
        entries.append({"accession": acc, "gtdb_leaf": ar_lookup[acc],
                        "tree": "ar53", "tokenized_path":
                        f"{TOKENIZED_BASE}/{acc}.npy"})

    random.shuffle(entries)

    # Optionally enrich with taxonomy
    if args.manifest and os.path.exists(args.manifest):
        print(f"Loading taxonomy from {args.manifest}...")
        tax_map = {}
        with open(args.manifest, newline="") as f:
            for row in csv.DictReader(f):
                tax_map[row["accession"]] = row
        for e in entries:
            row = tax_map.get(e["accession"], {})
            e["domain"] = row.get("domain", "")
            e["phylum"] = row.get("phylum", "")
            e["genus"]  = row.get("genus", "")

    accessions = [e["accession"] for e in entries]
    leaves_ordered = [e["gtdb_leaf"] for e in entries]

    # ── Parse trees (only once per tree) ──────────────────────────────────
    bac_indices = [i for i, e in enumerate(entries) if e["tree"] == "bac120"]
    ar_indices  = [i for i, e in enumerate(entries) if e["tree"] == "ar53"]

    # Full [N x N] distance matrix — cross-tree distances set to max
    N = len(entries)
    D_full = np.full((N, N), np.nan, dtype=np.float32)
    np.fill_diagonal(D_full, 0.0)

    if bac_indices:
        print(f"\nParsing bac120 tree ({len(bac_nwk)//1000}kB)...", flush=True)
        bac_tree = parse_tree(bac_nwk)
        bac_leaves = [leaves_ordered[i] for i in bac_indices]
        D_bac = compute_patristic_distances(bac_tree, bac_leaves, None)
        for ii, i in enumerate(bac_indices):
            for jj, j in enumerate(bac_indices):
                D_full[i, j] = D_bac[ii, jj]

    if ar_indices:
        print(f"\nParsing ar53 tree ({len(ar_nwk)//1000}kB)...", flush=True)
        ar_tree = parse_tree(ar_nwk)
        ar_leaves_list = [leaves_ordered[i] for i in ar_indices]
        D_ar = compute_patristic_distances(ar_tree, ar_leaves_list, None)
        for ii, i in enumerate(ar_indices):
            for jj, j in enumerate(ar_indices):
                D_full[i, j] = D_ar[ii, jj]

    # Cross-tree pairs: use a large sentinel value (max observed within-tree)
    within_max = np.nanmax(D_full[~np.isnan(D_full)])
    cross_dist = within_max * 1.5  # outside the ball
    D_full = np.where(np.isnan(D_full), cross_dist, D_full)

    # Normalize to [0, 1]
    d_max = D_full.max()
    if d_max > 0:
        D_norm = D_full / d_max
    else:
        D_norm = D_full.copy()

    print(f"\nPatristic distance stats:")
    print(f"  Raw   range: [{D_full[D_full>0].min():.6f}, {d_max:.6f}]")
    print(f"  Normalized:  [0, 1] by dividing by {d_max:.6f}")

    # ── Save outputs ──────────────────────────────────────────────────────
    manifest_path = out_dir / "telescope_manifest.csv"
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "gtdb_leaf", "tree",
                                                "tokenized_path",
                                                "domain", "phylum", "genus"])
        writer.writeheader()
        for e in entries:
            row = {k: e.get(k, "") for k in writer.fieldnames}
            writer.writerow(row)

    np.save(out_dir / "telescope_patristic.npy", D_norm)

    with open(out_dir / "telescope_patristic_raw.npy", "wb") as f:
        np.save(f, D_full)

    leaves_path = out_dir / "telescope_leaves.txt"
    with open(leaves_path, "w") as f:
        for acc in accessions:
            f.write(acc + "\n")

    print(f"\nSaved to {out_dir}/:")
    print(f"  telescope_manifest.csv      ({N} entries)")
    print(f"  telescope_patristic.npy     ({N}x{N}, normalized)")
    print(f"  telescope_patristic_raw.npy ({N}x{N}, raw branch lengths)")
    print(f"  telescope_leaves.txt        ({N} accessions)")
    print(f"\nNext: scp {out_dir}/ rohit@100.86.142.125:~/telescope/")
    print(f"Then: python3 fit_kappa_telescope.py --telescope-dir ~/telescope/")


if __name__ == "__main__":
    main()
