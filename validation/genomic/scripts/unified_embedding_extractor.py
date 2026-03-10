#!/usr/bin/env python3
"""
Robust Unified Hyperbolic Embedding Extractor
=============================================

Extracts hyperbolic coordinates from multiple BiosphereCodec models using
their specific config.json files to ensure perfect model instantiation.
This is the definitive method for extracting comparable embeddings.
"""

import os
import sys
import json
import numpy as np
import torch
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Add parent directory to path to allow BiosphereCodec import
sys.path.append(str(Path(__file__).parent))
try:
    from BiosphereCodec import BiosphereCodec
except ImportError:
    logger.error("Could not import BiosphereCodec. Make sure it's in the same directory.")
    sys.exit(1)


def load_model_from_config(model_dir: Path, checkpoint_file: Path) -> Optional[torch.nn.Module]:
    """Loads a model using its dedicated config.json file."""
    config_file = model_dir / 'config.json'
    if not config_file.is_file():
        logger.error(f"  ❌ Config file not found in {model_dir}")
        return None
    
    logger.info(f"  Loading model using config: {config_file}")
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # Extract relevant model parameters from config
    model_params = config.get('model_params', {})
    
    # Use the correct constructor argument name `vocab`
    vocab_size = model_params.get('vocab_size', 5444) 
    d_model = model_params.get('d_model', 256)
    n_layers = model_params.get('n_layers', 4)
    latent_dim = model_params.get('latent_dim', 256)
    max_len = model_params.get('max_len', 8192)

    logger.info(f"  Instantiating model with: vocab={vocab_size}, d_model={d_model}, n_layers={n_layers}, latent_dim={latent_dim}, max_len={max_len}")
    
    try:
        model = BiosphereCodec(
            vocab=vocab_size, 
            d_model=d_model, 
            n_layers=n_layers, 
            latent_dim=latent_dim,
            max_len=max_len
        )
        
        # Load the state dict
        checkpoint = torch.load(checkpoint_file, map_location='cpu')
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()
        logger.info("  ✅ Model instantiated and weights loaded successfully.")
        return model

    except Exception as e:
        logger.error(f"  ❌ Failed to build or load model from {model_dir}: {e}", exc_info=True)
        return None

def extract_embeddings(model: torch.nn.Module, checkpoint_file: Path) -> Tuple[Optional[np.ndarray], Optional[List[str]]]:
    """Extracts embeddings and genome IDs from a loaded model."""
    try:
        # Correct way to get the full embedding table
        embeddings = model.encoder.embed.weight.detach().cpu().numpy()
        checkpoint = torch.load(checkpoint_file, map_location='cpu')
        # The 'genome_ids' in the checkpoint likely correspond to the batch, not the full vocab
        # For the full embedding table, we generate IDs based on vocab index
        genome_ids = [f"token_{i}" for i in range(embeddings.shape[0])]
        logger.info(f"  ✅ Embeddings extracted. Shape: {embeddings.shape}")
        return embeddings, genome_ids
    except Exception as e:
        logger.error(f"  ❌ Failed to extract embeddings: {e}", exc_info=True)
        return None, None


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Robust Unified Hyperbolic Embedding Extractor")
    parser.add_argument("--base_dir", type=str, default="./data/genomes", help="Base directory of model runs.")
    parser.add_argument("--checkpoint_step", type=int, default=7000, help="Common checkpoint step to use.")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the output .npy files.")
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    models_to_process = {
        'original': 'biosphere_run_final',
        'seed42': 'biosphere_run_seed_42',
        'seed137': 'biosphere_run_seed_137',
        'seed2024': 'biosphere_run_seed_2024',
        'seed888': 'biosphere_run_seed_888'
    }

    all_genome_ids = None
    all_embeddings_extracted = True

    for name, model_folder in models_to_process.items():
        logger.info("-" * 50)
        logger.info(f"Processing model: {name}")
        model_dir = base_dir / model_folder
        checkpoint_file = model_dir / f"checkpoint_step_{args.checkpoint_step}.pt"

        if not checkpoint_file.exists():
            logger.error(f"Checkpoint file not found, skipping: {checkpoint_file}")
            all_embeddings_extracted = False
            continue

        model = load_model_from_config(model_dir, checkpoint_file)
        if model:
            embeddings, genome_ids = extract_embeddings(model, checkpoint_file)
            if embeddings is not None:
                output_filename = output_dir / f"emb_{name}_7k.npy"
                np.save(output_filename, embeddings)
                logger.info(f"  💾 Saved embeddings to {output_filename}")

                if all_genome_ids is None and genome_ids is not None:
                    all_genome_ids = genome_ids
                    ids_filename = output_dir / "genome_ids_7k.npy"
                    np.save(ids_filename, np.array(all_genome_ids, dtype=object))
                    logger.info(f"  💾 Saved genome IDs to {ids_filename}")
            else:
                all_embeddings_extracted = False
        else:
            all_embeddings_extracted = False

    logger.info("=" * 50)
    if all_embeddings_extracted:
        logger.info("🎉🎉🎉 Unified extraction complete for all models at step 7000.")
    else:
        logger.warning("⚠️ Some models could not be processed. Please check the logs.")


if __name__ == "__main__":
    main() 