#!/usr/bin/env python3
import re
import hashlib
import statistics as stats
from collections import Counter, defaultdict
from typing import List

import numpy as np
import torch

from biosphere_codec_model import BiosphereCodec
from rna_tokenizer import RNATokenizer, load_fasta_sequences


SUB_RE = re.compile(r"subtype\s*([A-Z][0-9]?)", re.IGNORECASE)
ALT_RE = re.compile(r"HIV[-\s]?1[^>\n]*\b([A-Z][0-9]?)\b")


def extract_subtype(header: str, seq: str) -> str:
    m = SUB_RE.search(header)
    if m:
        return m.group(1).upper()
    m = ALT_RE.search(header)
    if m:
        return m.group(1).upper()
    # fallback: hash bucket
    return f"h{int(hashlib.md5(seq.encode()).hexdigest(), 16) % 10}"


def run_sweep(
    data_path: str = "./data/viruses/HIV1_multi.fasta",
    max_seq: int = 12000,
    max_len: int = 4096,
    batch: int = 4,
    kappas: List[float] = [1.25, 1.30, 1.32, 1.35, 1.40, 1.45, 1.55],
    seeds: List[int] = [123, 321, 777],
    min_per_label: int = 40,
):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tok = RNATokenizer("nucleotide")

    raw = load_fasta_sequences(data_path, max_seq)
    seqs: List[str] = []
    labels_raw: List[str] = []
    for h, s in raw:
        if len(s) >= 1000:
            s = s[:max_len]
            seqs.append(s)
            labels_raw.append(extract_subtype(h, s))

    counts = Counter(labels_raw)
    keep = {l for l, c in counts.items() if c >= min_per_label}
    seqs = [s for s, l in zip(seqs, labels_raw) if l in keep]
    labels = [l for l in labels_raw if l in keep]
    if len(set(labels)) < 2:
        print("Insufficient subtype diversity after filtering; falling back to hash5.")
        labels = [int(hashlib.md5(s.encode()).hexdigest(), 16) % 5 for s in seqs]
    else:
        # map to ints for InfoNCE
        uniq = sorted(set(labels))
        lab_to_int = {l: i for i, l in enumerate(uniq)}
        labels = [lab_to_int[l] for l in labels]

    print("Device:", device)
    print("Loaded", len(seqs), "sequences; labels:", dict(Counter(labels)))

    def balanced_batches(rng: np.random.Generator, max_steps: int = 150):
        lbl_to_idx = defaultdict(list)
        for i, l in enumerate(labels):
            lbl_to_idx[l].append(i)
        valid = [l for l, idx in lbl_to_idx.items() if len(idx) >= 2]
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
            gb = [torch.tensor([0, L // 2, L - 1], dtype=torch.long) for _ in range(batch)]
            yield X.to(device), y.to(device), gb
            steps += 1

    results = []
    for kappa in kappas:
        hex_vals = []
        for sd in seeds:
            rng = np.random.default_rng(sd)
            model = BiosphereCodec(
                vocab=tok.vocab_size_actual,
                d_model=384,
                n_layers=4,
                latent_dim=96,
                mask_id=tok.mask_token_id,
                hex_weight=0.8,
                dist_weight=0.0,
            ).to(device)
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
                    if step < 2:
                        print(f"seed {sd} kappa {kappa:.2f} labels:", y.cpu().numpy())
            hex_avg = hex_sum / max(nb, 1)
            hex_vals.append(hex_avg)
        mean = stats.mean(hex_vals)
        sd = stats.pstdev(hex_vals) if len(hex_vals) > 1 else 0.0
        results.append((kappa, mean, sd))
        print(f"kappa={kappa:.2f}  HEX(mean±sd)={mean:.4f}±{sd:.4f}")

    print("\n=== Ranked by HEX (lower is better) ===")
    for kappa, mean, sd in sorted(results, key=lambda x: x[1]):
        print(f"kappa={kappa:.2f}  HEX={mean:.4f}±{sd:.4f}")


if __name__ == "__main__":
    run_sweep()







