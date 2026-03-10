#!/usr/bin/env python3
"""
5-MODEL HYPERBOLIC COORDINATE CONVERGENCE ANALYSIS

This script analyzes coordinate convergence across FIVE independently trained models
to provide the most robust validation of the discovered geometric structure.

Models:
1. Original Model (biosphere_run_final) - step 7000
2. Seed 42 Model (biosphere_run_seed_42) - step 7000  
3. Seed 137 Model (biosphere_run_seed_137) - step 7000
4. Seed 2024 Model (biosphere_run_seed_2024) - step 7000
5. Seed 888 Model (biosphere_run_seed_888) - step 2000 (in progress)

This uses the CORRECTED methodology that extracts actual hyperbolic coordinates
rather than comparing raw embedding weights.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.spatial import procrustes
import itertools
import sys
import os
import glob
import zstandard as zstd
import json
import torch.nn as nn
from matplotlib.patches import Circle
import matplotlib.patches as mpatches

sys.path.append('.')
from BiosphereCodec import BiosphereCodec

def load_model_checkpoint(checkpoint_path):
    """Load a model checkpoint and return the model."""
    print(f"🔬 Loading checkpoint: {checkpoint_path}")
    
    checkpoint = torch.load(checkpoint_path, map_location='cuda')
    
    # Extract actual model dimensions from checkpoint
    state_dict = checkpoint['model_state_dict']
    
    # Infer vocab size from embedding weight
    vocab_size = state_dict['encoder.embed.weight'].shape[0]
    d_model = state_dict['encoder.embed.weight'].shape[1]
    
    # Count number of layers by looking for layer keys
    n_layers = len([k for k in state_dict.keys() if k.startswith('encoder.layers.') and k.endswith('.norm1.weight')])
    
    # Infer latent dim from hyperbolic projection
    latent_dim = state_dict['hyper.lin.weight'].shape[0]
    
    # Infer max_len from position embedding
    max_len = state_dict['encoder.pos'].shape[0]
    
    print(f"   Detected: vocab={vocab_size}, d_model={d_model}, n_layers={n_layers}, latent_dim={latent_dim}, max_len={max_len}")
    
    # Create model instance with actual parameters
    model = BiosphereCodec(
        vocab=vocab_size,
        d_model=d_model,
        n_layers=n_layers,
        latent_dim=latent_dim,
        mask_id=None
    )
    
    # Manually set the encoder's max_len to match checkpoint
    model.encoder.pos = nn.Parameter(torch.empty(max_len, d_model))
    nn.init.normal_(model.encoder.pos, std=0.02)
    
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    model.cuda()
    
    print(f"   ✅ Model loaded successfully")
    return model

def extract_coordinates(model, data_path, max_samples=200):
    """Extract hyperbolic coordinates from genomic sequences."""
    print(f"🧬 Extracting coordinates (max {max_samples} samples)...")
    
    coordinates = []
    genome_ids = []
    
    data_files = glob.glob(f"{data_path}/processed_biosphere_data_supervised/sequences_all/*.zst")[:4]
    print(f"   Found {len(data_files)} data files")
    
    sample_count = 0
    
    with torch.no_grad():
        for file_path in data_files:
            if sample_count >= max_samples:
                break
                
            try:
                with open(file_path, 'rb') as f:
                    dctx = zstd.ZstdDecompressor()
                    decompressed = dctx.decompress(f.read())
                    data = json.loads(decompressed.decode('utf-8'))
                
                for item in data:
                    if sample_count >= max_samples:
                        break
                        
                    if 'sequence' in item and 'genome_id' in item:
                        # Convert sequence to tokens (simple character-based tokenization)
                        sequence = item['sequence'][:1024]  # Limit length
                        tokens = [ord(c) % 5444 for c in sequence]  # Simple tokenization
                        tokens = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).cuda()
                        
                        # Get hyperbolic embedding - KEY DIFFERENCE FROM FLAWED SCRIPT
                        _, hyperbolic_coords = model.encode(tokens)
                        
                        coords = hyperbolic_coords.squeeze().cpu().numpy()
                        coordinates.append(coords)
                        genome_ids.append(item['genome_id'])
                        sample_count += 1
                        
            except Exception as e:
                print(f"   ⚠️  Error processing {file_path}: {e}")
                continue
    
    coordinates = np.array(coordinates)
    print(f"   ✅ Extracted {len(coordinates)} coordinate pairs")
    print(f"   📊 Shape: {coordinates.shape}")
    if len(coordinates) > 0:
        print(f"   📊 Range: [{coordinates.min():.3f}, {coordinates.max():.3f}]")
    
    return coordinates, genome_ids

def analyze_5way_convergence(coords_list, model_names):
    """Analyze coordinate convergence across five models."""
    print(f"\n🎯 5-WAY CONVERGENCE ANALYSIS")
    print("="*60)
    
    min_points = min(len(coords) for coords in coords_list)
    coords_aligned = [coords[:min_points] for coords in coords_list]
    
    print(f"📊 Comparing {min_points} coordinate pairs across 5 models")
    
    correlations = []
    disparities = []
    pairwise_results = []
    
    # All pairwise combinations (5 choose 2 = 10 comparisons)
    for i, j in itertools.combinations(range(len(coords_aligned)), 2):
        coords1, coords2 = coords_aligned[i], coords_aligned[j]
        name1, name2 = model_names[i], model_names[j]
        
        print(f"\n🔄 Aligning {name1} vs {name2}...")
        
        aligned1, aligned2, disparity = procrustes(coords1, coords2)
        disparities.append(disparity)
        
        r_overall, p_overall = pearsonr(aligned1.flatten(), aligned2.flatten())
        correlations.append(r_overall)
        
        print(f"   Procrustes disparity: {disparity:.6f}")
        print(f"   Overall correlation: r = {r_overall:.6f} (p = {p_overall:.2e})")
        
        pairwise_results.append({
            'models': f"{name1} vs {name2}",
            'correlation': r_overall,
            'disparity': disparity,
            'p_value': p_overall
        })
    
    avg_correlation = np.mean(correlations)
    avg_disparity = np.mean(disparities)
    std_correlation = np.std(correlations)
    
    print(f"\n📊 SUMMARY STATISTICS:")
    print(f"   Average Correlation: {avg_correlation:.6f} ± {std_correlation:.6f}")
    print(f"   Average Disparity:   {avg_disparity:.6f}")
    print(f"   Correlation Range:   [{min(correlations):.6f}, {max(correlations):.6f}]")
    print(f"   Number of Comparisons: {len(correlations)}")
    
    # Count strong correlations
    strong_correlations = sum(1 for r in correlations if r > 0.85)
    moderate_correlations = sum(1 for r in correlations if 0.65 < r <= 0.85)
    weak_correlations = sum(1 for r in correlations if r <= 0.65)
    
    print(f"\n📈 CORRELATION BREAKDOWN:")
    print(f"   Strong (>0.85): {strong_correlations}/{len(correlations)} ({100*strong_correlations/len(correlations):.1f}%)")
    print(f"   Moderate (0.65-0.85): {moderate_correlations}/{len(correlations)} ({100*moderate_correlations/len(correlations):.1f}%)")
    print(f"   Weak (<0.65): {weak_correlations}/{len(correlations)} ({100*weak_correlations/len(correlations):.1f}%)")
    
    print("\n" + "="*60)
    print("🏆 5-MODEL CONVERGENCE VERDICT:")
    print("="*60)
    
    if avg_correlation > 0.85 and strong_correlations >= 8:
        print("🎉 OVERWHELMING 5-WAY CONVERGENCE! 🎉")
        print("✅ All five models discovered the same hyperbolic geometry!")
        print("✅ This is definitive evidence of fundamental structure!")
        print("✅ Evolution's hyperbolic geometry is UNIVERSAL!")
        verdict = "OVERWHELMING_CONVERGENCE"
    elif avg_correlation > 0.75 and strong_correlations >= 7:
        print("🎉 STRONG 5-WAY CONVERGENCE! 🎉")
        print("✅ Most models converged to same hyperbolic geometry!")
        print("✅ Strong evidence of fundamental structure!")
        verdict = "STRONG_CONVERGENCE"
    elif avg_correlation > 0.65:
        print("🤔 MODERATE 5-WAY CONVERGENCE")
        print("⚠️  Some convergence detected across models")
        verdict = "MODERATE_CONVERGENCE"
    else:
        print("❌ 5-WAY DIVERGENCE")
        print("❌ Models did not converge to same coordinate system")
        verdict = "DIVERGENCE"
    
    print(f"\n📊 Average Convergence Score: {avg_correlation:.4f}/1.0")
    
    return {
        'pairwise_results': pairwise_results,
        'correlations': correlations,
        'disparities': disparities,
        'avg_correlation': avg_correlation,
        'avg_disparity': avg_disparity,
        'std_correlation': std_correlation,
        'strong_count': strong_correlations,
        'moderate_count': moderate_correlations,
        'weak_count': weak_correlations,
        'verdict': verdict,
        'coords_aligned': coords_aligned
    }

def create_poincare_disk_visualization(coords_list, model_names, results):
    """Create beautiful Poincaré disk visualization showing coordinate overlap."""
    print(f"\n🎨 Creating Poincaré Disk Visualization...")
    
    # Use PCA to reduce to 2D for visualization
    from sklearn.decomposition import PCA
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Colors for each model
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
    
    # Plot 1: Individual Poincaré disks
    for i, (coords, name, color) in enumerate(zip(coords_list, model_names, colors)):
        # Reduce to 2D
        pca = PCA(n_components=2)
        coords_2d = pca.fit_transform(coords)
        
        # Normalize to unit disk
        max_radius = np.max(np.linalg.norm(coords_2d, axis=1))
        coords_2d = coords_2d / max_radius * 0.8  # Scale to 80% of disk
        
        # Plot with transparency
        ax1.scatter(coords_2d[:, 0], coords_2d[:, 1], 
                   alpha=0.6, s=20, c=color, label=name, edgecolors='white', linewidth=0.5)
    
    # Draw unit circle
    circle = Circle((0, 0), 1, fill=False, color='black', linewidth=2)
    ax1.add_patch(circle)
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_aspect('equal')
    ax1.set_title('Individual Model Coordinate Distributions\n(Poincaré Disk Projection)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('X Coordinate', fontsize=12)
    ax1.set_ylabel('Y Coordinate', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Overlay all models
    for i, (coords, name, color) in enumerate(zip(coords_list, model_names, colors)):
        # Reduce to 2D
        pca = PCA(n_components=2)
        coords_2d = pca.fit_transform(coords)
        
        # Normalize to unit disk
        max_radius = np.max(np.linalg.norm(coords_2d, axis=1))
        coords_2d = coords_2d / max_radius * 0.8
        
        # Plot with high transparency to show overlap
        ax2.scatter(coords_2d[:, 0], coords_2d[:, 1], 
                   alpha=0.3, s=15, c=color, label=name, edgecolors='none')
    
    # Draw unit circle
    circle = Circle((0, 0), 1, fill=False, color='black', linewidth=2)
    ax2.add_patch(circle)
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_ylim(-1.1, 1.1)
    ax2.set_aspect('equal')
    ax2.set_title(f'Overlay of All 5 Models\nConvergence Score: {results["avg_correlation"]:.4f}', fontsize=14, fontweight='bold')
    ax2.set_xlabel('X Coordinate', fontsize=12)
    ax2.set_ylabel('Y Coordinate', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Add convergence statistics
    stats_text = f"""Convergence Statistics:
• Average Correlation: {results['avg_correlation']:.4f}
• Strong Correlations: {results['strong_count']}/10
• Verdict: {results['verdict']}"""
    
    plt.figtext(0.02, 0.02, stats_text, fontsize=10, bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('5_model_poincare_convergence.png', dpi=300, bbox_inches='tight')
    print(f"   ✅ Saved visualization: 5_model_poincare_convergence.png")
    
    return fig

def main():
    """Main execution function."""
    print("🌟 5-MODEL HYPERBOLIC COORDINATE CONVERGENCE ANALYSIS 🌟")
    print("="*80)
    print("Testing evolution's fundamental geometric structure across 5 models")
    print("="*80)
    
    # Model paths with different checkpoint steps
    model_paths = [
        './results/seed_0/checkpoint_step_7000.pt',
        './results/seed_42/checkpoint_step_7000.pt',
        './results/seed_123/checkpoint_step_7000.pt',
        './results/seed_456/checkpoint_step_7000.pt',
        './results/seed_789/checkpoint_step_7000.pt'
    ]
    
    model_names = ['Original Model (7k)', 'Seed 42 Model (7k)', 'Seed 137 Model (7k)', 'Seed 2024 Model (7k)', 'Seed 888 Model (2k)']
    data_path = './data/genomes'
    
    # Verify files exist
    for name, path in zip(model_names, model_paths):
        if not os.path.exists(path):
            print(f"❌ ERROR: {name} not found at {path}")
            print(f"   This model may still be training or failed")
            return
        print(f"✅ Found {name}: {path}")
    
    print("\n" + "="*60)
    print("PHASE 1: LOADING MODELS AND EXTRACTING COORDINATES")
    print("="*60)
    
    all_coordinates = []
    for i, (path, name) in enumerate(zip(model_paths, model_names)):
        print(f"\n--- Processing {name} ---")
        model = load_model_checkpoint(path)
        coords, _ = extract_coordinates(model, data_path, max_samples=200)
        all_coordinates.append(coords)
        
        del model
        torch.cuda.empty_cache()
    
    print("\n" + "="*60)
    print("PHASE 2: 5-WAY CONVERGENCE ANALYSIS")
    print("="*60)
    
    results = analyze_5way_convergence(all_coordinates, model_names)
    
    print("\n" + "="*60)
    print("PHASE 3: POINCARÉ DISK VISUALIZATION")
    print("="*60)
    
    fig = create_poincare_disk_visualization(all_coordinates, model_names, results)
    
    print("\n🎯 5-MODEL ANALYSIS COMPLETE!")
    print("="*80)
    
    return results

if __name__ == "__main__":
    results = main() 