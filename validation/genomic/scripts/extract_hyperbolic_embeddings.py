#!/usr/bin/env python3
"""
Extract hyperbolic (Poincaré-ball) embeddings from all 5 models for the same set of genomes.
Saves embeddings and genome_ids for alignment and visualization.
"""
import torch
import numpy as np
import glob
import json
import zstandard as zstd
import sys
import os

sys.path.append('/zfs_raid/SentryBio/5k_test_genomes')
from BiosphereCodec import BiosphereCodec

def load_model_checkpoint(checkpoint_path):
    checkpoint = torch.load(checkpoint_path, map_location='cuda')
    state_dict = checkpoint['model_state_dict']
    vocab_size = state_dict['encoder.embed.weight'].shape[0]
    d_model = state_dict['encoder.embed.weight'].shape[1]
    n_layers = len([k for k in state_dict.keys() if k.startswith('encoder.layers.') and k.endswith('.norm1.weight')])
    latent_dim = state_dict['hyper.lin.weight'].shape[0]
    max_len = state_dict['encoder.pos'].shape[0]
    model = BiosphereCodec(
        vocab=vocab_size,
        d_model=d_model,
        n_layers=n_layers,
        latent_dim=latent_dim,
        mask_id=None
    )
    model.encoder.pos = torch.nn.Parameter(torch.empty(max_len, d_model))
    torch.nn.init.normal_(model.encoder.pos, std=0.02)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    model.cuda()
    return model

def extract_embeddings(model, data_path, genome_ids_ref=None, max_samples=200):
    coordinates = []
    genome_ids = []
    data_files = sorted(glob.glob(f"{data_path}/processed_biosphere_data_supervised/sequences_all/*.zst"))
    sample_count = 0
    with torch.no_grad():
        for file_path in data_files:
            if sample_count >= max_samples:
                break
            with open(file_path, 'rb') as f:
                dctx = zstd.ZstdDecompressor()
                decompressed = dctx.decompress(f.read())
                data = json.loads(decompressed.decode('utf-8'))
            for item in data:
                if sample_count >= max_samples:
                    break
                if 'sequence' in item and 'genome_id' in item:
                    if genome_ids_ref is not None and item['genome_id'] not in genome_ids_ref:
                        continue
                    sequence = item['sequence'][:1024]
                    tokens = [ord(c) % 5444 for c in sequence]
                    tokens = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).cuda()
                    _, hyperbolic_coords = model.encode(tokens)
                    coords = hyperbolic_coords.squeeze().cpu().numpy()
                    coordinates.append(coords)
                    genome_ids.append(item['genome_id'])
                    sample_count += 1
    return np.array(coordinates), genome_ids

def main():
    data_path = '/zfs_raid/SentryBio/5k_test_genomes'
    model_ckpts = [
        ('emb_original.npy', '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_final/checkpoint_step_7000.pt'),
        ('emb_seed42.npy', '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_42/checkpoint_step_7000.pt'),
        ('emb_seed137.npy', '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_137/checkpoint_step_7000.pt'),
        ('emb_seed2024.npy', '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_2024/checkpoint_step_7000.pt'),
        ('emb_seed888.npy', '/zfs_raid/SentryBio/5k_test_genomes/biosphere_run_seed_888/checkpoint_step_7000.pt'),
    ]
    # Extract reference genome_ids from the first model
    print(f"Extracting reference embeddings and genome_ids from {model_ckpts[0][1]}")
    model = load_model_checkpoint(model_ckpts[0][1])
    coords, genome_ids = extract_embeddings(model, data_path, genome_ids_ref=None, max_samples=200)
    np.save('emb_original.npy', coords)
    np.save('genome_ids.npy', np.array(genome_ids))
    print(f"Saved emb_original.npy and genome_ids.npy ({len(genome_ids)} samples)")
    # For other models, use the same genome_ids (order must match)
    for fname, ckpt in model_ckpts[1:]:
        print(f"Extracting embeddings from {ckpt}")
        model = load_model_checkpoint(ckpt)
        coords, _ = extract_embeddings(model, data_path, genome_ids_ref=genome_ids, max_samples=len(genome_ids))
        np.save(fname, coords)
        print(f"Saved {fname} ({coords.shape[0]} samples)")

if __name__ == "__main__":
    main() 