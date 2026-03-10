#!/usr/bin/env python3
"""
Null Controls for κ Validation
================================
1. Shuffled labels: Random permutation → κ preference should flatten
2. Single label: All sequences same label → HEX ≈ 0, κ uninformative
3. Random sequences: Synthetic data → no phylogenetic signal
"""

import hashlib
import statistics as stats
from collections import Counter, defaultdict
from typing import List
import numpy as np
import torch
from biosphere_codec_model import BiosphereCodec
from rna_tokenizer import RNATokenizer, load_fasta_sequences


def run_null_control_sweep(
    data_path: str,
    control_type: str,
    name: str,
    max_seq: int = 5000,
    max_len: int = 4096,
    batch: int = 4,
    kappas: List[float] = [1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60],
    seeds: List[int] = [13, 17, 23],
):
    """
    Run κ sweep with null control labeling.
    
    control_type:
        - 'shuffled': Random permutation of authentic labels
        - 'single': All sequences assigned label 0
        - 'random_hash': Random hash labels (no phylogenetic relationship)
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tok = RNATokenizer("nucleotide")

    # Load sequences
    raw = load_fasta_sequences(data_path, max_seq)
    seqs: List[str] = []
    authentic_labels: List[int] = []
    
    for h, s in raw:
        if len(s) >= 1000:
            s = s[:max_len]
            seqs.append(s)
            # Create authentic hash5 labels for reference
            authentic_labels.append(int(hashlib.md5(s.encode()).hexdigest(), 16) % 5)
    
    print(f"\n{'='*60}")
    print(f"NULL CONTROL: {control_type.upper()}")
    print(f"Dataset: {name}")
    print(f"{'='*60}")
    print(f"Device: {device}")
    print(f"Loaded {len(seqs)} sequences")
    
    # Apply null control labeling
    if control_type == "shuffled":
        print(f"Authentic label distribution: {dict(Counter(authentic_labels))}")
        # Shuffle labels randomly (breaks phylogenetic association)
        rng_shuffle = np.random.default_rng(42)
        labels = authentic_labels.copy()
        rng_shuffle.shuffle(labels)
        print(f"Shuffled label distribution: {dict(Counter(labels))}")
        print("→ Expect: κ preference should flatten (no real hierarchy)")
    
    elif control_type == "single":
        # All sequences same label → no contrastive pairs
        labels = [0] * len(seqs)
        print(f"Single label: all sequences → label 0")
        print("→ Expect: HEX ≈ 0 (no pairs), κ uninformative")
    
    elif control_type == "random_hash":
        # Random hash labels (uniform distribution, no phylogeny)
        rng_hash = np.random.default_rng(123)
        labels = [int(rng_hash.integers(0, 5)) for _ in seqs]
        print(f"Random hash distribution: {dict(Counter(labels))}")
        print("→ Expect: κ preference should be weak/flat")
    
    else:
        raise ValueError(f"Unknown control_type: {control_type}")
    
    def balanced_batches(rng: np.random.Generator, max_steps: int = 150):
        """Generate balanced batches (or fail gracefully for single-label)."""
        lbl_to_idx = defaultdict(list)
        for i, l in enumerate(labels):
            lbl_to_idx[l].append(i)
        
        valid = [l for l, idx in lbl_to_idx.items() if len(idx) >= 2]
        
        if len(valid) < 2:
            print(f"⚠️  Only {len(valid)} valid labels → cannot form pairs, HEX will be 0")
            # Fallback: yield random batches (will have HEX=0)
            for step in range(max_steps):
                picks = rng.choice(len(seqs), size=batch, replace=False)
                ids = [torch.tensor(tok.encode(seqs[k], add_special_tokens=True)[:max_len]) for k in picks]
                L = max(x.size(0) for x in ids)
                X = torch.zeros(batch, L, dtype=torch.long)
                for i, x in enumerate(ids):
                    X[i, : x.size(0)] = x
                y = torch.tensor([labels[k] for k in picks], dtype=torch.long)
                gb = [torch.tensor([0, L // 2, L - 1], dtype=torch.long) for _ in range(batch)]
                yield X.to(device), y.to(device), gb
            return
        
        # Normal balanced batching
        steps = 0
        while steps < max_steps and len(valid) >= 2:
            l1, l2 = rng.choice(valid, 2, replace=False)
            a = lbl_to_idx[l1][:]
            b = lbl_to_idx[l2][:]
            rng.shuffle(a)
            rng.shuffle(b)
            if len(a) < 2 or len(b) < 2:
                continue
            picks = [a[0], a[1], b[0], b[1]]
            ids = [torch.tensor(tok.encode(seqs[k], add_special_tokens=True)[:max_len]) for k in picks]
            L = max(x.size(0) for x in ids)
            X = torch.zeros(batch, L, dtype=torch.long)
            for i, x in enumerate(ids):
                X[i, : x.size(0)] = x
            y = torch.tensor([labels[k] for k in picks], dtype=torch.long)
            binc = np.bincount(y.cpu().numpy())
            if not (binc >= 2).any():
                continue
            gb = [torch.tensor([0, L // 2, L - 1], dtype=torch.long) for _ in range(batch)]
            yield X.to(device), y.to(device), gb
            steps += 1
    
    # Run sweep
    results = []
    for kappa in kappas:
        hex_vals = []
        for sd in seeds:
            rng = np.random.default_rng(sd)
            model = BiosphereCodec(
                vocab=tok.vocab_size_actual,
                d_model=256,
                n_layers=3,
                latent_dim=64,
                mask_id=tok.mask_token_id,
                hex_weight=0.8,
                dist_weight=0.0,
            ).to(device)
            
            # Pin κ
            with torch.no_grad():
                if hasattr(model.hyper, "manifold") and hasattr(model.hyper.manifold, "c"):
                    model.hyper.manifold.c.data.fill_(kappa)
                elif hasattr(model.hyper, "c"):
                    model.hyper.c.data.fill_(kappa)
            
            for n, p in model.named_parameters():
                if n.endswith(".c"):
                    p.requires_grad = False
            
            model.loss_fn.temp = 0.05
            
            hex_sum = 0.0
            nb = 0
            with torch.no_grad():
                for step, (X, y, gb) in enumerate(balanced_batches(rng, 150)):
                    _, logs = model(X, gene_idx=gb, tax_ids=y)
                    hex_sum += float(logs.get("hex", 0.0))
                    nb += 1
            
            hex_avg = hex_sum / max(nb, 1)
            hex_vals.append(hex_avg)
        
        mean = stats.mean(hex_vals)
        sd_val = stats.pstdev(hex_vals) if len(hex_vals) > 1 else 0.0
        results.append((kappa, mean, sd_val))
        print(f"κ={kappa:.2f}  HEX(mean±sd)={mean:.4f}±{sd_val:.4f}")
    
    print("\n=== Ranked by HEX (lower is better) ===")
    sorted_results = sorted(results, key=lambda x: x[1])
    for kappa, mean, sd_val in sorted_results:
        print(f"κ={kappa:.2f}  HEX={mean:.4f}±{sd_val:.4f}")
    
    # Compute κ preference strength (range of HEX values)
    hex_values = [r[1] for r in results]
    hex_range = max(hex_values) - min(hex_values)
    print(f"\nHEX range: {hex_range:.4f}")
    
    if control_type == "single":
        print(f"Average HEX: {stats.mean(hex_values):.4f} (expect ≈0 for single-label)")
    
    return sorted_results, hex_range


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--data_path", type=str, required=True)
    p.add_argument("--control_type", type=str, required=True, 
                   choices=["shuffled", "single", "random_hash"])
    p.add_argument("--name", type=str, required=True)
    p.add_argument("--max_sequences", type=int, default=5000)
    args = p.parse_args()
    
    run_null_control_sweep(
        data_path=args.data_path,
        control_type=args.control_type,
        name=args.name,
        max_seq=args.max_sequences,
    )


if __name__ == "__main__":
    main()


