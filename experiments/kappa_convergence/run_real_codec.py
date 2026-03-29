#!/usr/bin/env python3
"""
Run the REAL BiosphereCodec.py from June 30 2025 — UNCHANGED.
Only this file is new: a minimal data loader + training loop.
"""

import argparse, json, logging, math, os, random, sys, time
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from torch.optim import AdamW

# Import the REAL BiosphereCodec — no modifications
sys.path.insert(0, os.path.dirname(__file__))
CODEC_DIR = str(Path(__file__).resolve().parent.parent.parent.parent
                / "Biosphere_codec" / "external" / "nexus_py")
sys.path.insert(0, CODEC_DIR)
from BiosphereCodec import BiosphereCodec


class NpyGenomeDataset(Dataset):
    """Load .npy tokenized genomes with genus labels for hex_loss."""

    def __init__(self, manifest_path, max_len=8192, max_genomes=None,
                 vocab_size=4096, min_genus_count=2):
        import pandas as pd
        df = pd.read_csv(manifest_path, low_memory=False)
        df = df[df['tokenized_path'].notna() & df['genus'].notna()]
        df = df[df['tokenized_path'].apply(lambda p: os.path.exists(str(p)))]

        gc = df['genus'].value_counts()
        valid = gc[gc >= min_genus_count].index
        df = df[df['genus'].isin(valid)]

        if max_genomes and len(df) > max_genomes:
            df = df.sample(max_genomes, random_state=42)

        self.paths = df['tokenized_path'].tolist()
        genera = df['genus'].tolist()
        self.genus_to_id = {g: i for i, g in enumerate(sorted(set(genera)))}
        self.labels = [self.genus_to_id[g] for g in genera]
        self.max_len = max_len
        self.vocab_size = vocab_size
        print(f"Dataset: {len(self.paths)} genomes, {len(self.genus_to_id)} genera")

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        tokens = np.load(self.paths[idx]).astype(np.int64)
        tokens = np.clip(tokens, 0, self.vocab_size - 1)
        if len(tokens) > self.max_len:
            start = random.randint(0, len(tokens) - self.max_len)
            tokens = tokens[start:start + self.max_len]
        return torch.from_numpy(tokens), self.labels[idx]


def collate_fn(batch):
    tokens_list, labels = zip(*batch)
    max_len = max(t.shape[0] for t in tokens_list)
    padded = torch.zeros(len(tokens_list), max_len, dtype=torch.long)
    for i, t in enumerate(tokens_list):
        padded[i, :t.shape[0]] = t
    tax_ids = torch.tensor(labels, dtype=torch.long)
    return padded, tax_ids


def train(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s',
                        handlers=[logging.FileHandler(outdir / "training.log"),
                                  logging.StreamHandler()])
    log = logging.getLogger("real_codec")

    # === THE REAL BiosphereCodec — geoopt intact, nothing changed ===
    model = BiosphereCodec(
        vocab=args.vocab_size,
        d_model=args.d_model,
        n_layers=args.n_layers,
        max_len=args.max_len,
        latent_dim=args.latent_dim,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"Model: {n_params/1e6:.1f}M params, d={args.d_model}, "
             f"L={args.n_layers}, latent={args.latent_dim}")
    log.info(f"Initial c = {model.hyper.c.item():.6f}")
    log.info(f"c.requires_grad = {model.hyper.c.requires_grad}")
    log.info(f"Using REAL BiosphereCodec.py from {CODEC_DIR}")

    dataset = NpyGenomeDataset(
        args.manifest, max_len=args.max_len,
        max_genomes=args.max_genomes, vocab_size=args.vocab_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True,
                        collate_fn=collate_fn, num_workers=0, drop_last=True)

    optimizer = AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)

    log.info(f"Training: {args.steps} steps, batch={args.batch_size}, device={device}")

    model.train()
    step = 0
    kappa_history = []
    loader_iter = iter(loader)

    while step < args.steps:
        try:
            tokens, tax_ids = next(loader_iter)
        except StopIteration:
            loader_iter = iter(loader)
            tokens, tax_ids = next(loader_iter)

        tokens = tokens.to(device)
        tax_ids = tax_ids.to(device)

        # THE KEY: pass tax_ids to activate Poincaré InfoNCE (hex_loss)
        loss, logs = model(tokens, tax_ids=tax_ids)

        optimizer.zero_grad()
        loss.backward()

        # Read gradient BEFORE clipping/stepping
        c_grad = model.hyper.c.grad.item() if model.hyper.c.grad is not None else 0.0

        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        step += 1

        kappa = model.hyper.c.item()
        kappa_history.append({'step': step, 'kappa': kappa, 'c_grad': c_grad, **logs})

        if step % args.log_every == 0:
            log.info(f"Step {step:>5d} | c={kappa:.6f} | dc={c_grad:+.4e} | "
                     f"MLM={logs['mlm']:.4f} hex={logs['hex']:.4f} "
                     f"dist={logs['dist']:.4f}")

        if step % 1000 == 0:
            torch.save({
                'step': step, 'model_state_dict': model.state_dict(),
                'kappa': kappa, 'kappa_history': kappa_history[-200:],
            }, outdir / f"checkpoint_{step}.pt")

    final_kappa = model.hyper.c.item()
    torch.save({'step': step, 'model_state_dict': model.state_dict(),
                'kappa': final_kappa}, outdir / "final.pt")
    with open(outdir / "kappa_history.json", 'w') as f:
        json.dump(kappa_history, f, indent=2)

    log.info(f"COMPLETE: Final c = {final_kappa:.6f}")
    return final_kappa


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--output-dir", default="./real_codec_run")
    p.add_argument("--vocab-size", type=int, default=4096)
    p.add_argument("--d-model", type=int, default=256)
    p.add_argument("--n-layers", type=int, default=4)
    p.add_argument("--latent-dim", type=int, default=256)
    p.add_argument("--max-len", type=int, default=8192)
    p.add_argument("--max-genomes", type=int, default=5000)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--lr", type=float, default=2e-4)
    p.add_argument("--steps", type=int, default=7000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log-every", type=int, default=10)
    p.add_argument("--sweep", action="store_true")
    args = p.parse_args()

    if args.sweep:
        seeds = [0, 42, 123, 456, 789]
        results = []
        base = args.output_dir
        for seed in seeds:
            args.seed = seed
            args.output_dir = f"{base}/seed_{seed}"
            k = train(args)
            results.append({'seed': seed, 'kappa': k})
            with open(f"{base}/sweep.json", 'w') as f:
                json.dump(results, f, indent=2)
        kappas = [r['kappa'] for r in results]
        mean_k, std_k = np.mean(kappas), np.std(kappas)
        print(f"\nκ = {mean_k:.6f} ± {std_k:.6f} (CV={std_k/mean_k*100:.1f}%)")
    else:
        train(args)


if __name__ == "__main__":
    main()
