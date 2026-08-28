#!/usr/bin/env python3
"""
E11: Pure 2D Coordinate System (Layer IIa instrument)

WHAT THIS TESTS: Whether a polar H^2 chart — kappa frozen, no classification
heads, no ODE — is seed-stable up to O(2). That is reproducibility within
the imposed model (H^2, frozen kappa, taxonomy quartets, radial target),
not a host-class selection test and not a derivation of kappa from the
genetic code. Radius-as-depth and angle-as-divergence are modeling choices.

ARCHITECTURE:
  Token embedding -> mean pool -> MLP -> exp_map -> H^2(kappa frozen)
  ~35K parameters. No classification heads. No ODE.
  kappa is a design constant. InfoNCE temperature is degenerate with
  curvature, so kappa cannot be discovered by gradient descent. The
  reference run uses 5/4; that number is not a theorem.

LOSS (three terms):
  1. Domain separation: organisms from different domains should be far apart
     in Poincare distance. Simple contrastive margin loss.
  2. Quartet consistency: d(a,b) + d(c,d) < d(a,c) + d(b,d)
     Quartets derived from taxonomy (genus < family < order < class < phylum < domain).
     This is the four-point classifier used as a training signal, not a
     curvature calibrator.
  3. Radial ordering: larger genome -> deeper in ball.
     Genome size is an independent depth proxy. It is not accumulated
     information and not a clock (E6: CCS radial axis stays advisory).

ANCHOR SELECTION:
  200-300 organisms from the full manifest, balanced across:
  - All three domains
  - Major phyla within each domain
  - Full range of genome sizes

SEED STABILITY TEST:
  Train N seeds. Procrustes-align the resulting coordinate sets.
  If the polar chart is determined by the data inside this imposed model,
  residual after alignment is small compared with ball radius; the leftover
  gauge is global O(2) (rotation and reflection) unless orientation is
  fixed. This does not compare host classes and does not certify a filled
  atlas of life.

USAGE:
  python train.py \\
    --manifest /path/to/manifest.csv \\
    --output-dir /path/to/E11_results \\
    --device cuda \\
    --n-seeds 5
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

# ── Constants ──────────────────────────────────────────────────────────────

KAPPA = 1.25          # 5/4 — frozen design constant, not a learned parameter.
VOCAB_SIZE = 4096
FEAT_DIM = 8
MAX_SEQ_LEN = 8192    # tokens per genome window
LATENT_DIM = 2        # the whole point

TAXONOMY_RANKS = ["genus", "family", "order", "class", "phylum", "domain"]

# Landmark organisms for WGS84 cross-seed comparison
LANDMARK_ORGANISMS = [
    {"name": "Homo sapiens",                  "domain": "Eukaryota", "genus": "Homo",                "species_match": "Homo sapiens"},
    {"name": "Saccharomyces cerevisiae",      "domain": "Eukaryota", "genus": "Saccharomyces",      "species_match": "Saccharomyces cerevisiae"},
    {"name": "Arabidopsis thaliana",          "domain": "Eukaryota", "genus": "Arabidopsis",         "species_match": "Arabidopsis thaliana"},
    {"name": "Escherichia coli",              "domain": "Bacteria",  "genus": "Escherichia",         "species_match": "Escherichia coli"},
    {"name": "Bacillus subtilis",             "domain": "Bacteria",  "genus": "Bacillus",            "species_match": "Bacillus subtilis"},
    {"name": "Methanocaldococcus jannaschii", "domain": "Archaea",   "genus": "Methanocaldococcus",  "species_match": "Methanocaldococcus jannaschii"},
    {"name": "Halobacterium salinarum",       "domain": "Archaea",   "genus": "Halobacterium",       "species_match": "Halobacterium salinarum"},
]


# ── Poincare Ball Operations (curvature = -KAPPA) ─────────────────────────

def poincare_distance(x: torch.Tensor, y: torch.Tensor, c: float = KAPPA) -> torch.Tensor:
    """Geodesic distance in the Poincare ball with curvature -c."""
    sqrt_c = math.sqrt(c)
    x_sq = (x * x).sum(-1, keepdim=True).clamp(max=1.0 / c - 1e-5)
    y_sq = (y * y).sum(-1, keepdim=True).clamp(max=1.0 / c - 1e-5)
    xy_diff_sq = ((x - y) ** 2).sum(-1, keepdim=True)
    num = 2.0 * c * xy_diff_sq
    denom = (1.0 - c * x_sq) * (1.0 - c * y_sq)
    arg = 1.0 + num / denom.clamp(min=1e-10)
    return ((1.0 / sqrt_c) * torch.acosh(arg.clamp(min=1.0 + 1e-7))).squeeze(-1)


def exp_map_zero(v: torch.Tensor, c: float = KAPPA) -> torch.Tensor:
    """Exponential map at the origin of the Poincare ball."""
    sqrt_c = math.sqrt(c)
    v_norm = v.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    return torch.tanh(sqrt_c * v_norm / 2.0) * v / (sqrt_c * v_norm)


def project_to_ball(x: torch.Tensor, c: float = KAPPA, eps: float = 1e-5) -> torch.Tensor:
    """Project onto open ball {x : c||x||^2 < 1}."""
    max_norm = (1.0 / math.sqrt(c)) - eps
    norms = x.norm(dim=-1, keepdim=True).clamp(min=1e-8)
    clamped = (max_norm / norms).clamp(max=1.0)
    return x * clamped


# ── Model ──────────────────────────────────────────────────────────────────

class PoincareEncoder2D(nn.Module):
    """Minimal encoder: tokens -> 2D Poincare ball coordinates.

    Architecture: Embedding(4096,8) -> Conv1d(8,32,k=5) -> GELU
                  -> adaptive mean pool -> Linear(32,64) -> GELU
                  -> Linear(64,64) -> GELU -> Linear(64,2) -> exp_map

    The conv layer captures local k-mer motifs that distinguish clades.
    Total: ~38K parameters.
    """

    CONV_DIM = 32

    def __init__(self, vocab_size: int = VOCAB_SIZE, feat_dim: int = FEAT_DIM,
                 hidden: int = 64):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, feat_dim)
        self.conv = nn.Conv1d(feat_dim, self.CONV_DIM, kernel_size=5, padding=2)
        self.encoder = nn.Sequential(
            nn.Linear(self.CONV_DIM, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, LATENT_DIM),
        )

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        """tokens: (B, L) long -> coordinates: (B, 2) on Poincare ball."""
        x = self.embed(tokens)          # (B, L, feat_dim)
        x = x.transpose(1, 2)          # (B, feat_dim, L) for conv
        x = F.gelu(self.conv(x))        # (B, CONV_DIM, L)
        x = x.mean(dim=2)              # (B, CONV_DIM) — global mean pool
        v = self.encoder(x)             # (B, 2) — tangent vector at origin
        z = exp_map_zero(v)             # (B, 2) — on the ball
        return project_to_ball(z)


# ── Dataset ────────────────────────────────────────────────────────────────

class AnchorDataset(Dataset):
    """Dataset for anchor organisms with random-window augmentation."""

    def __init__(self, anchors: List[Dict], max_len: int = MAX_SEQ_LEN,
                 augment: bool = True):
        self.anchors = anchors
        self.max_len = max_len
        self.augment = augment

    def __len__(self) -> int:
        return len(self.anchors)

    def __getitem__(self, idx: int) -> Dict:
        a = self.anchors[idx]
        path = a["tokenized_path"]

        if path.endswith(".npy"):
            tokens = torch.from_numpy(np.load(path)).long()
        else:
            tokens = torch.load(path, map_location="cpu", weights_only=False)
            if isinstance(tokens, dict):
                tokens = tokens.get("input_ids", tokens.get("tokens"))
            tokens = tokens.long().squeeze()

        # Augmentation: random contiguous window
        if self.augment and len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start : start + self.max_len]
        elif len(tokens) > self.max_len:
            tokens = tokens[: self.max_len]

        return {
            "tokens": tokens,
            "idx": idx,
            "domain": a.get("domain", "Unknown"),
            "phylum": a.get("phylum", "Unknown"),
            "class": a.get("class", "Unknown"),
            "order": a.get("order", "Unknown"),
            "family": a.get("family", "Unknown"),
            "genus": a.get("genus", "Unknown"),
            "genome_size": float(a.get("genome_size", 0) or 0),
        }


def collate_anchors(batch: List[Dict]) -> Dict:
    """Pad variable-length token sequences and collate metadata."""
    max_len = max(b["tokens"].shape[0] for b in batch)
    padded = torch.zeros(len(batch), max_len, dtype=torch.long)
    for i, b in enumerate(batch):
        L = b["tokens"].shape[0]
        padded[i, :L] = b["tokens"]

    return {
        "tokens": padded,
        "idx": torch.tensor([b["idx"] for b in batch]),
        "domains": [b["domain"] for b in batch],
        "phyla": [b["phylum"] for b in batch],
        "classes": [b["class"] for b in batch],
        "orders": [b["order"] for b in batch],
        "families": [b["family"] for b in batch],
        "genera": [b["genus"] for b in batch],
        "genome_sizes": torch.tensor(
            [b["genome_size"] for b in batch], dtype=torch.float32
        ),
    }


# ── Anchor Selection ──────────────────────────────────────────────────────

def select_anchors(manifest_path: str, n_target: int = 250,
                   seed: int = 42) -> List[Dict]:
    """Select balanced anchor organisms from manifest.

    Strategy: group by domain -> phylum -> family, then take evenly-spaced
    samples within each family sorted by genome size to maximize diversity.
    """
    rng = random.Random(seed)

    all_rows = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            tp = row.get("tokenized_path", "")
            if not tp or not os.path.exists(tp):
                continue
            if row.get("domain") in ("Viruses", ""):
                continue
            all_rows.append(row)

    # Group by domain -> phylum -> family
    tree: Dict[str, Dict[str, Dict[str, List]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for row in all_rows:
        d = row.get("domain", "Unknown")
        p = row.get("phylum", "Unknown")
        f = row.get("family", "Unknown")
        tree[d][p][f].append(row)

    domains = sorted(tree.keys())
    n_domains = len(domains)

    # Phase 1: proportional allocation with per-domain minimums
    domain_budget = {}
    for d in domains:
        n_available = sum(
            len(fam) for phy in tree[d].values() for fam in phy.values()
        )
        # At least 30 per domain, at most proportional share
        domain_budget[d] = max(30, int(n_target * n_available / len(all_rows)))

    anchors = []
    for domain in domains:
        budget = domain_budget[domain]
        phyla = sorted(tree[domain].keys())
        n_per_phylum = max(2, budget // max(len(phyla), 1))

        for phylum in phyla:
            families = sorted(tree[domain][phylum].keys())
            n_per_family = max(1, n_per_phylum // max(len(families), 1))

            for family in families:
                candidates = tree[domain][phylum][family]
                candidates.sort(key=lambda r: float(r.get("genome_size", 0) or 0))

                if len(candidates) <= n_per_family:
                    anchors.extend(candidates)
                else:
                    indices = np.linspace(
                        0, len(candidates) - 1, n_per_family, dtype=int
                    )
                    anchors.extend([candidates[i] for i in indices])

    # Phase 2: trim if grossly over target
    if len(anchors) > n_target * 2:
        rng.shuffle(anchors)
        # Keep at least 30 per domain
        trimmed = []
        by_domain = defaultdict(list)
        for a in anchors:
            by_domain[a.get("domain", "Unknown")].append(a)
        for d in domains:
            pool = by_domain[d]
            take = max(30, int(n_target * len(pool) / len(anchors)))
            trimmed.extend(pool[:take])
        anchors = trimmed

    return anchors


# ── Taxonomy Distance ─────────────────────────────────────────────────────

def taxonomy_distance(a: Dict, b: Dict) -> int:
    """Distance between two samples based on shared taxonomic ranks.

    Returns 0 (same genus) through 6 (different domain).
    Lower = more closely related.
    """
    for level, rank in enumerate(TAXONOMY_RANKS):
        va = a.get(rank, "Unknown")
        vb = b.get(rank, "Unknown")
        if va == vb and va != "Unknown" and va != "":
            return level
    return len(TAXONOMY_RANKS)


def batch_taxonomy_distance(batch: Dict, i: int, j: int) -> int:
    """Taxonomy distance between two indices within a batch."""
    for level, rank_key in enumerate(
        ["genera", "families", "orders", "classes", "phyla", "domains"]
    ):
        va = batch[rank_key][i]
        vb = batch[rank_key][j]
        if va == vb and va not in ("Unknown", ""):
            return level
    return 6


# ── Loss Functions ─────────────────────────────────────────────────────────

def quartet_loss(coords: torch.Tensor, batch: Dict,
                 n_quartets: int = 200, margin: float = 0.5) -> torch.Tensor:
    """Quartet consistency loss derived from taxonomy.

    For each quartet (a,b,c,d), find the topology consistent with taxonomy
    and enforce the four-point condition:
        d(a,b) + d(c,d) < d(a,c) + d(b,d) + margin
    """
    B = coords.shape[0]
    if B < 4:
        return torch.tensor(0.0, device=coords.device)

    indices = list(range(B))
    n_quartets = min(n_quartets, B * 3)
    loss = torch.tensor(0.0, device=coords.device)
    count = 0

    for _ in range(n_quartets):
        q = random.sample(indices, 4)

        # Three possible splits: (01|23), (02|13), (03|12)
        splits = [
            (q[0], q[1], q[2], q[3]),
            (q[0], q[2], q[1], q[3]),
            (q[0], q[3], q[1], q[2]),
        ]

        # Correct split: minimize sum of within-pair taxonomy distances
        best = min(
            splits,
            key=lambda s: (
                batch_taxonomy_distance(batch, s[0], s[1])
                + batch_taxonomy_distance(batch, s[2], s[3])
            ),
        )
        i, j, k, l = best

        # Four-point condition: d(i,j) + d(k,l) should be smallest
        d_ij = poincare_distance(coords[i : i + 1], coords[j : j + 1])
        d_kl = poincare_distance(coords[k : k + 1], coords[l : l + 1])
        d_ik = poincare_distance(coords[i : i + 1], coords[k : k + 1])
        d_jl = poincare_distance(coords[j : j + 1], coords[l : l + 1])

        violation = (d_ij + d_kl) - (d_ik + d_jl) + margin
        loss = loss + F.relu(violation).squeeze()
        count += 1

    return loss / max(count, 1)


def domain_separation_loss(coords: torch.Tensor, domains: List[str],
                           n_pairs: int = 300, margin: float = 0.5) -> torch.Tensor:
    """Different domains should be far apart; same domain should be close.

    Simple contrastive loss on Poincare distances:
    - Same domain: d(i,j) should be small (no penalty if < margin_same)
    - Different domain: d(i,j) should be large (penalty if < margin_cross)
    """
    B = coords.shape[0]
    if B < 2:
        return torch.tensor(0.0, device=coords.device)

    indices = list(range(B))
    loss = torch.tensor(0.0, device=coords.device)
    count = 0

    for _ in range(n_pairs):
        i, j = random.sample(indices, 2)
        di, dj = domains[i], domains[j]
        if di in ("Unknown", "") or dj in ("Unknown", ""):
            continue

        d = poincare_distance(coords[i:i+1], coords[j:j+1]).squeeze()

        if di != dj:
            # Cross-domain: want distance > margin
            loss = loss + F.relu(margin - d)
        else:
            # Same domain: want distance < margin * 0.3
            loss = loss + F.relu(d - margin * 0.3)
        count += 1

    return loss / max(count, 1)


def radial_ordering_loss(coords: torch.Tensor, genome_sizes: torch.Tensor,
                         n_pairs: int = 300, margin: float = 0.02) -> torch.Tensor:
    """Larger genome -> deeper in ball (higher r).

    Samples pairs where one genome is at least 1.5x larger than the other
    and penalizes if the smaller genome has a larger radius.
    """
    r = coords.norm(dim=-1)
    B = r.shape[0]
    if B < 2:
        return torch.tensor(0.0, device=coords.device)

    # Filter to valid genome sizes
    valid = (genome_sizes > 0).nonzero(as_tuple=True)[0].tolist()
    if len(valid) < 2:
        return torch.tensor(0.0, device=coords.device)

    n_pairs = min(n_pairs, len(valid) * 3)
    loss = torch.tensor(0.0, device=coords.device)
    count = 0

    for _ in range(n_pairs):
        i, j = random.sample(valid, 2)
        si, sj = genome_sizes[i].item(), genome_sizes[j].item()

        if si > sj * 1.5:
            # i should be deeper (higher r)
            loss = loss + F.relu(r[j] - r[i] + margin)
            count += 1
        elif sj > si * 1.5:
            loss = loss + F.relu(r[i] - r[j] + margin)
            count += 1

    return loss / max(count, 1)


# ── Training ───────────────────────────────────────────────────────────────

def train_one_seed(
    anchors: List[Dict],
    seed: int,
    n_epochs: int,
    lr: float,
    batch_size: int,
    device: str,
    log,
) -> Tuple[nn.Module, np.ndarray]:
    """Train one model from a given seed. Return model and final coordinates."""
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    model = PoincareEncoder2D().to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    dataset = AnchorDataset(anchors, augment=True)
    loader = DataLoader(
        dataset,
        batch_size=min(len(anchors), batch_size),
        shuffle=True,
        collate_fn=collate_anchors,
        num_workers=0,
    )

    n_params = sum(p.numel() for p in model.parameters())
    log(f"  Seed {seed}: {n_params:,} parameters")

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        epoch_lq = 0.0
        epoch_ld = 0.0
        epoch_lr = 0.0
        n_batches = 0

        for batch in loader:
            tokens = batch["tokens"].to(device)
            coords = model(tokens)

            lq = quartet_loss(coords, batch)
            ld = domain_separation_loss(coords, batch["domains"])
            lr_loss = radial_ordering_loss(
                coords, batch["genome_sizes"].to(device)
            )

            loss = lq + ld + 0.3 * lr_loss

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_loss += loss.item()
            epoch_lq += lq.item()
            epoch_ld += ld.item()
            epoch_lr += lr_loss.item()
            n_batches += 1

        scheduler.step()

        if (epoch + 1) % 50 == 0 or epoch == 0 or epoch == n_epochs - 1:
            avg = epoch_loss / max(n_batches, 1)
            aq = epoch_lq / max(n_batches, 1)
            ad = epoch_ld / max(n_batches, 1)
            ar = epoch_lr / max(n_batches, 1)
            log(
                f"    Epoch {epoch + 1:4d}: loss={avg:.4f}  "
                f"(quartet={aq:.4f}  domain={ad:.4f}  radial={ar:.4f})"
            )

    # Extract final coordinates (no augmentation, deterministic)
    model.eval()
    eval_dataset = AnchorDataset(anchors, augment=False)
    eval_loader = DataLoader(
        eval_dataset,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=collate_anchors,
        num_workers=0,
    )

    all_coords = []
    with torch.no_grad():
        for batch in eval_loader:
            tokens = batch["tokens"].to(device)
            coords = model(tokens)
            all_coords.append(coords.cpu())

    return model, torch.cat(all_coords, dim=0).numpy()


# ── Evaluation ─────────────────────────────────────────────────────────────

def evaluate_coordinates(
    coords: np.ndarray, anchors: List[Dict], log
) -> Dict:
    """Evaluate quality of the learned 2D coordinates."""
    r = np.sqrt((coords ** 2).sum(axis=1))
    theta = np.arctan2(coords[:, 1], coords[:, 0])

    # 1. Quartet consistency rate
    n_correct = 0
    n_tested = 0
    indices = list(range(len(anchors)))

    rng = random.Random(123)
    for _ in range(min(2000, len(anchors) * 5)):
        q = rng.sample(indices, 4)
        splits = [
            (q[0], q[1], q[2], q[3]),
            (q[0], q[2], q[1], q[3]),
            (q[0], q[3], q[1], q[2]),
        ]

        # Correct split from taxonomy
        correct = min(
            splits,
            key=lambda s: (
                taxonomy_distance(anchors[s[0]], anchors[s[1]])
                + taxonomy_distance(anchors[s[2]], anchors[s[3]])
            ),
        )

        # Geometric split: which split minimizes d(i,j) + d(k,l)?
        def split_cost(s):
            d1 = np.sqrt(((coords[s[0]] - coords[s[1]]) ** 2).sum())
            d2 = np.sqrt(((coords[s[2]] - coords[s[3]]) ** 2).sum())
            return d1 + d2

        geometric = min(splits, key=split_cost)

        if geometric == correct:
            n_correct += 1
        n_tested += 1

    quartet_rate = n_correct / max(n_tested, 1)
    log(f"  Quartet consistency: {n_correct}/{n_tested} ({quartet_rate:.1%})")

    # 2. Domain centroid separation
    domain_centroids = {}
    for a, c in zip(anchors, coords):
        d = a.get("domain", "Unknown")
        if d not in domain_centroids:
            domain_centroids[d] = []
        domain_centroids[d].append(c)

    for d in sorted(domain_centroids):
        pts = np.array(domain_centroids[d])
        centroid = pts.mean(axis=0)
        spread = np.sqrt(((pts - centroid) ** 2).sum(axis=1).mean())
        log(f"  {d:12s}: centroid=({centroid[0]:+.4f}, {centroid[1]:+.4f})  "
            f"spread={spread:.4f}  n={len(pts)}")

    domains = list(domain_centroids.keys())
    domain_dists = {}
    for i in range(len(domains)):
        ci = np.array(domain_centroids[domains[i]]).mean(axis=0)
        for j in range(i + 1, len(domains)):
            cj = np.array(domain_centroids[domains[j]]).mean(axis=0)
            dist = np.sqrt(((ci - cj) ** 2).sum())
            pair = f"{domains[i][:1]}-{domains[j][:1]}"
            domain_dists[pair] = dist
            log(f"  {domains[i]}-{domains[j]} distance: {dist:.4f}")

    # 3. Radial correlation with genome size
    genome_sizes = np.array([float(a.get("genome_size", 0) or 0) for a in anchors])
    valid = genome_sizes > 0
    if valid.sum() > 10:
        from scipy.stats import spearmanr

        rho, pval = spearmanr(r[valid], genome_sizes[valid])
        log(f"  Radial-genome_size Spearman: rho={rho:.4f}  p={pval:.2e}")
    else:
        rho, pval = 0.0, 1.0
        log(f"  Radial-genome_size: insufficient data")

    # 4. r distribution
    log(f"  r range: [{r.min():.4f}, {r.max():.4f}]  mean={r.mean():.4f} +/- {r.std():.4f}")

    return {
        "quartet_consistency": quartet_rate,
        "domain_distances": domain_dists,
        "radial_genome_size_rho": float(rho),
        "radial_genome_size_p": float(pval),
        "r_min": float(r.min()),
        "r_max": float(r.max()),
        "r_mean": float(r.mean()),
    }


# ── Seed Stability Test (WGS84) ───────────────────────────────────────────

def procrustes_align(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, float]:
    """Align Y to X via Procrustes (rotation + reflection). Return aligned Y
    and mean residual distance."""
    mx, my = X.mean(0), Y.mean(0)
    Xc, Yc = X - mx, Y - my

    # Normalize scale
    sx = np.sqrt((Xc ** 2).sum())
    sy = np.sqrt((Yc ** 2).sum())
    Xc /= max(sx, 1e-10)
    Yc /= max(sy, 1e-10)

    M = Xc.T @ Yc
    U, S, Vt = np.linalg.svd(M)
    R = Vt.T @ U.T

    # Handle reflection
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T

    Y_aligned = (Yc @ R) * sx + mx
    residual = np.sqrt(((X - Y_aligned) ** 2).sum(axis=1).mean())

    return Y_aligned, residual


def seed_stability_test(
    all_coords: List[np.ndarray], log
) -> Dict:
    """Compare coordinate sets across seeds using Procrustes alignment."""
    n_seeds = len(all_coords)
    if n_seeds < 2:
        log("  Only 1 seed — skipping stability test.")
        return {"verdict": "SINGLE_SEED"}

    residuals = np.zeros((n_seeds, n_seeds))

    for i in range(n_seeds):
        for j in range(i + 1, n_seeds):
            _, res = procrustes_align(all_coords[i], all_coords[j])
            residuals[i, j] = res
            residuals[j, i] = res

    triu = residuals[np.triu_indices(n_seeds, k=1)]
    mean_res = triu.mean()
    max_res = triu.max()
    min_res = triu.min()

    log(f"  Pairwise Procrustes residuals ({n_seeds} seeds):")
    log(f"    Mean:  {mean_res:.6f}")
    log(f"    Min:   {min_res:.6f}")
    log(f"    Max:   {max_res:.6f}")

    # Full matrix
    log(f"  Residual matrix:")
    seeds_labels = [f"S{i}" for i in range(n_seeds)]
    header = "       " + "  ".join(f"{s:>8s}" for s in seeds_labels)
    log(header)
    for i in range(n_seeds):
        row = f"  {seeds_labels[i]:>4s}  "
        row += "  ".join(
            f"{residuals[i, j]:8.5f}" if i != j else "    ---  "
            for j in range(n_seeds)
        )
        log(row)

    # WGS84 verdict
    # Ball radius = 1/sqrt(kappa) = 1/sqrt(1.25) ~ 0.894
    ball_radius = 1.0 / math.sqrt(KAPPA)
    threshold = 0.05 * ball_radius  # 5% of ball radius

    if mean_res < threshold:
        verdict = "COORDINATE_STABLE"
        log(f"\n  VERDICT: {verdict}")
        log(f"  Mean residual {mean_res:.6f} < threshold {threshold:.6f}")
        log(f"  Coordinates are reproducible across random seeds.")
        log(f"  This is a datum, not a model artifact.")
    elif mean_res < threshold * 2:
        verdict = "COORDINATE_MARGINAL"
        log(f"\n  VERDICT: {verdict}")
        log(f"  Mean residual {mean_res:.6f} near threshold {threshold:.6f}")
        log(f"  Coordinates are approximately stable. May improve with more epochs.")
    else:
        verdict = "COORDINATE_UNSTABLE"
        log(f"\n  VERDICT: {verdict}")
        log(f"  Mean residual {mean_res:.6f} > threshold {threshold:.6f}")
        log(f"  Coordinates vary with initialization.")
        log(f"  Consider: more anchors, more epochs, or adjusted loss weights.")

    return {
        "mean_residual": float(mean_res),
        "max_residual": float(max_res),
        "min_residual": float(min_res),
        "threshold": float(threshold),
        "ball_radius": float(ball_radius),
        "verdict": verdict,
    }


# ── Landmark Organism Evaluation ───────────────────────────────────────────

def find_landmark_organisms(manifest_path: str) -> List[Dict]:
    """Find one representative accession per landmark organism in manifest."""
    landmarks = []
    with open(manifest_path) as f:
        rows = list(csv.DictReader(f))

    for lm in LANDMARK_ORGANISMS:
        match = None
        for row in rows:
            species = row.get("species", "")
            genus = row.get("genus", "")
            tp = row.get("tokenized_path", "")
            if not tp or not os.path.exists(tp):
                continue
            if genus == lm["genus"] and lm["species_match"] in species:
                match = row
                break
        if match:
            landmarks.append({**match, "landmark_name": lm["name"]})
        else:
            landmarks.append({"landmark_name": lm["name"], "accession": "NOT_FOUND"})

    return landmarks


def landmark_cross_seed_table(
    landmarks: List[Dict], models: List[nn.Module], device: str, log
) -> Dict:
    """Run each landmark through all seed models and report (r, theta) table."""
    n_seeds = len(models)
    found = [lm for lm in landmarks if lm.get("accession") != "NOT_FOUND"]

    if not found:
        log("  No landmark organisms found in manifest.")
        return {}

    log(f"  Found {len(found)}/{len(landmarks)} landmark organisms")
    log()

    # Header
    seed_headers = "  ".join(f"{'Seed ' + str(i):>16s}" for i in range(n_seeds))
    log(f"  {'Organism':<32s}  {'Domain':<10s}  {seed_headers}  {'r std':>8s}  {'theta std':>10s}")
    log(f"  {'-'*32}  {'-'*10}  " + "  ".join(["-" * 16] * n_seeds) + f"  {'-'*8}  {'-'*10}")

    table_data = []
    for lm in found:
        path = lm["tokenized_path"]
        if path.endswith(".npy"):
            tokens = torch.from_numpy(np.load(path)).long()
        else:
            tokens = torch.load(path, map_location="cpu", weights_only=False)
            if isinstance(tokens, dict):
                tokens = tokens.get("input_ids", tokens.get("tokens"))
            tokens = tokens.long().squeeze()

        # Truncate to MAX_SEQ_LEN (no augmentation — deterministic)
        if len(tokens) > MAX_SEQ_LEN:
            tokens = tokens[:MAX_SEQ_LEN]

        tokens_batch = tokens.unsqueeze(0).to(device)

        rs = []
        thetas = []
        coords_all = []
        for model in models:
            model.eval()
            with torch.no_grad():
                coord = model(tokens_batch).cpu().numpy()[0]
            r = float(np.sqrt(coord[0] ** 2 + coord[1] ** 2))
            theta = float(np.degrees(np.arctan2(coord[1], coord[0])))
            rs.append(r)
            thetas.append(theta)
            coords_all.append(coord)

        r_std = np.std(rs)
        theta_std = np.std(thetas)

        seed_strs = "  ".join(
            f"({rs[i]:.4f}, {thetas[i]:+6.1f}deg)" for i in range(n_seeds)
        )
        name = lm["landmark_name"]
        domain = lm.get("domain", "?")
        log(f"  {name:<32s}  {domain:<10s}  {seed_strs}  {r_std:8.5f}  {theta_std:+9.2f} deg")

        table_data.append({
            "name": name,
            "accession": lm.get("accession", ""),
            "domain": domain,
            "r_per_seed": rs,
            "theta_per_seed": thetas,
            "r_mean": float(np.mean(rs)),
            "r_std": float(r_std),
            "theta_mean": float(np.mean(thetas)),
            "theta_std": float(theta_std),
        })

    # Summary
    log()
    all_r_stds = [t["r_std"] for t in table_data]
    all_theta_stds = [t["theta_std"] for t in table_data]
    log(f"  Mean r std across landmarks:     {np.mean(all_r_stds):.5f}")
    log(f"  Mean theta std across landmarks: {np.mean(all_theta_stds):.2f} deg")
    log(f"  Max r std:                       {np.max(all_r_stds):.5f}")
    log(f"  Max theta std:                   {np.max(all_theta_stds):.2f} deg")

    return {"landmarks": table_data}


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="E11: Pure 2D Coordinate System"
    )
    parser.add_argument("--manifest", required=True,
                        help="Path to manifest CSV with tokenized_path column")
    parser.add_argument("--output-dir", default="E11_results")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--n-anchors", type=int, default=250)
    parser.add_argument("--n-seeds", type=int, default=5)
    parser.add_argument("--n-epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, "E11_full.log")
    log_file = open(log_path, "w")

    def log(msg=""):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}" if msg else ""
        print(line)
        log_file.write(line + "\n")
        log_file.flush()

    log("=" * 70)
    log("E11: PURE 2D COORDINATE SYSTEM")
    log("=" * 70)
    log(f"kappa = {KAPPA} (frozen — theorem, not parameter)")
    log(f"latent_dim = {LATENT_DIM}")
    log(f"Manifest: {args.manifest}")
    log(f"Target anchors: {args.n_anchors}")
    log(f"Seeds: {args.n_seeds}")
    log(f"Epochs per seed: {args.n_epochs}")
    log()

    # ── Select anchors ────────────────────────────────────────────────────
    log("--- Selecting anchor organisms ---")
    anchors = select_anchors(args.manifest, n_target=args.n_anchors)

    domain_counts = defaultdict(int)
    phylum_counts = defaultdict(int)
    family_set = set()
    for a in anchors:
        domain_counts[a.get("domain", "Unknown")] += 1
        phylum_counts[a.get("phylum", "Unknown")] += 1
        family_set.add(a.get("family", "Unknown"))

    log(f"  Selected {len(anchors)} anchors:")
    for d, c in sorted(domain_counts.items()):
        log(f"    {d}: {c}")
    log(f"    {len(phylum_counts)} phyla, {len(family_set)} families")

    # Genome size range
    gsizes = [float(a.get("genome_size", 0) or 0) for a in anchors]
    valid_gsizes = [g for g in gsizes if g > 0]
    if valid_gsizes:
        log(
            f"    Genome sizes: {min(valid_gsizes)/1e6:.1f} - "
            f"{max(valid_gsizes)/1e6:.1f} Mb "
            f"(median {np.median(valid_gsizes)/1e6:.1f} Mb)"
        )
    log()

    # ── Train multiple seeds ──────────────────────────────────────────────
    all_coords = []
    all_models = []

    for seed_idx in range(args.n_seeds):
        seed = 42 + seed_idx * 17
        log(f"--- Training seed {seed_idx + 1}/{args.n_seeds} (seed={seed}) ---")
        model, coords = train_one_seed(
            anchors, seed, args.n_epochs, args.lr, args.batch_size,
            args.device, log,
        )
        all_coords.append(coords)
        all_models.append(model)

        r = np.sqrt((coords ** 2).sum(axis=1))
        log(f"    r range: [{r.min():.4f}, {r.max():.4f}], mean={r.mean():.4f}")
        log()

    # ── Evaluate reference coordinates (seed 0) ──────────────────────────
    log("=" * 70)
    log("COORDINATE EVALUATION (reference seed)")
    log("=" * 70)
    eval_results = evaluate_coordinates(all_coords[0], anchors, log)
    log()

    # ── Seed stability test ───────────────────────────────────────────────
    log("=" * 70)
    log("SEED STABILITY TEST (WGS84)")
    log("=" * 70)
    stability = seed_stability_test(all_coords, log)
    log()

    # ── Landmark organism cross-seed comparison ──────────────────────────
    log("=" * 70)
    log("LANDMARK ORGANISMS — (r, theta) ACROSS SEEDS")
    log("=" * 70)
    landmarks = find_landmark_organisms(args.manifest)
    landmark_results = landmark_cross_seed_table(
        landmarks, all_models, args.device, log
    )
    log()

    # ── Save results ──────────────────────────────────────────────────────
    coord_records = []
    ref_coords = all_coords[0]
    for i, a in enumerate(anchors):
        coord_records.append({
            "accession": a.get("accession", ""),
            "domain": a.get("domain", ""),
            "phylum": a.get("phylum", ""),
            "family": a.get("family", ""),
            "genus": a.get("genus", ""),
            "genome_size": float(a.get("genome_size", 0) or 0),
            "x": float(ref_coords[i, 0]),
            "y": float(ref_coords[i, 1]),
            "r": float(np.sqrt(ref_coords[i, 0] ** 2 + ref_coords[i, 1] ** 2)),
            "theta": float(np.arctan2(ref_coords[i, 1], ref_coords[i, 0])),
        })

    results = {
        "experiment": "E11",
        "kappa": KAPPA,
        "kappa_note": "frozen constant = 5/4 from state equation",
        "latent_dim": LATENT_DIM,
        "n_anchors": len(anchors),
        "n_seeds": args.n_seeds,
        "n_epochs": args.n_epochs,
        "learning_rate": args.lr,
        "domain_counts": dict(domain_counts),
        "evaluation": eval_results,
        "stability": stability,
        "landmarks": landmark_results.get("landmarks", []),
        "coordinates": coord_records,
    }

    # Convert numpy types for JSON serialization
    def to_json_safe(obj):
        if isinstance(obj, (np.floating, np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")

    out_path = os.path.join(args.output_dir, "E11_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=to_json_safe)
    log(f"Results saved to {out_path}")

    # Save reference model
    model_path = os.path.join(args.output_dir, "E11_model_seed0.pt")
    torch.save(all_models[0].state_dict(), model_path)
    log(f"Reference model saved to {model_path}")

    # Save all seed coordinates for analysis
    for i, coords in enumerate(all_coords):
        np.save(
            os.path.join(args.output_dir, f"coords_seed{i}.npy"), coords
        )

    log()
    log("=" * 70)
    log("E11 COMPLETE")
    log("=" * 70)

    log_file.close()


if __name__ == "__main__":
    main()
