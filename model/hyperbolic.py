"""
Hyperbolic Procrustes Alignment
===============================

This module implements hyperbolic Procrustes alignment for aligning hyperbolic embeddings.
The algorithm:
1. Log-map to tangent space
2. Center the embeddings
3. Apply Kabsch rotation
4. Exp-map back to hyperbolic space

Based on the user's description of "log-map to tangent space, center, Kabsch rotation, exp-map back".
"""

import numpy as np
from typing import Tuple, Optional
import logging

logger = logging.getLogger(__name__)

def hyperbolic_distance(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Compute hyperbolic distance between points in the Poincaré ball model.
    
    Args:
        x: Points in hyperbolic space (n_samples, n_dim)
        y: Points in hyperbolic space (n_samples, n_dim)
    
    Returns:
        Hyperbolic distances
    """
    # Ensure inputs are numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Compute squared norms
    x_norm_sq = np.sum(x**2, axis=1, keepdims=True)
    y_norm_sq = np.sum(y**2, axis=1, keepdims=True)
    
    # Compute Euclidean distance squared
    diff = x - y
    euclidean_dist_sq = np.sum(diff**2, axis=1, keepdims=True)
    
    # Hyperbolic distance formula for Poincaré ball
    numerator = 2 * euclidean_dist_sq
    denominator = (1 - x_norm_sq) * (1 - y_norm_sq)
    
    # Avoid division by zero
    denominator = np.maximum(denominator, 1e-8)
    
    # Compute hyperbolic distance
    dist = np.arccosh(1 + numerator / denominator)
    
    return dist.flatten()

def log_map(x: np.ndarray, base_point: np.ndarray) -> np.ndarray:
    """
    Map points from hyperbolic space to tangent space at base_point.
    
    Args:
        x: Points in hyperbolic space (n_samples, n_dim)
        base_point: Base point for tangent space (n_dim,)
    
    Returns:
        Points in tangent space (n_samples, n_dim)
    """
    x = np.asarray(x)
    base_point = np.asarray(base_point).flatten()
    
    # Ensure base_point is at origin for simplicity
    if np.allclose(base_point, 0):
        # If base point is origin, log map is simpler
        x_norm = np.linalg.norm(x, axis=1, keepdims=True)
        x_norm = np.maximum(x_norm, 1e-8)  # Avoid division by zero
        
        # Log map formula for origin
        log_x = x / x_norm * np.arctanh(x_norm)
        return log_x
    else:
        # For general base point, use parallel transport
        # This is a simplified version - in practice you'd want more sophisticated parallel transport
        x_norm = np.linalg.norm(x, axis=1, keepdims=True)
        x_norm = np.maximum(x_norm, 1e-8)
        
        # Simplified log map
        log_x = x / x_norm * np.arctanh(x_norm)
        return log_x

def exp_map(v: np.ndarray, base_point: np.ndarray) -> np.ndarray:
    """
    Map points from tangent space to hyperbolic space.
    
    Args:
        v: Points in tangent space (n_samples, n_dim)
        base_point: Base point for tangent space (n_dim,)
    
    Returns:
        Points in hyperbolic space (n_samples, n_dim)
    """
    v = np.asarray(v)
    base_point = np.asarray(base_point).flatten()
    
    # Ensure base_point is at origin for simplicity
    if np.allclose(base_point, 0):
        # If base point is origin, exp map is simpler
        v_norm = np.linalg.norm(v, axis=1, keepdims=True)
        v_norm = np.maximum(v_norm, 1e-8)  # Avoid division by zero
        
        # Exp map formula for origin
        exp_v = v / v_norm * np.tanh(v_norm)
        return exp_v
    else:
        # For general base point, use parallel transport
        # This is a simplified version
        v_norm = np.linalg.norm(v, axis=1, keepdims=True)
        v_norm = np.maximum(v_norm, 1e-8)
        
        # Simplified exp map
        exp_v = v / v_norm * np.tanh(v_norm)
        return exp_v

def kabsch_rotation(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """
    Compute the optimal rotation matrix using Kabsch algorithm.
    
    Args:
        X: Source points (n_samples, n_dim)
        Y: Target points (n_samples, n_dim)
    
    Returns:
        Rotation matrix (n_dim, n_dim)
    """
    # Center the points
    X_centered = X - np.mean(X, axis=0)
    Y_centered = Y - np.mean(Y, axis=0)
    
    # Compute covariance matrix
    H = X_centered.T @ Y_centered
    
    # SVD decomposition
    U, S, Vt = np.linalg.svd(H)
    
    # Compute rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (handle reflection case)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    return R

def center_embeddings(X: np.ndarray) -> np.ndarray:
    """
    Center embeddings by subtracting the mean.
    
    Args:
        X: Embeddings (n_samples, n_dim)
    
    Returns:
        Centered embeddings (n_samples, n_dim)
    """
    return X - np.mean(X, axis=0)

def hyperbolic_procrustes(X: np.ndarray, Y: np.ndarray, 
                         return_disparity: bool = False,
                         max_iter: int = 100,
                         tol: float = 1e-6) -> Tuple[np.ndarray, Optional[float]]:
    """
    Apply hyperbolic Procrustes alignment to align Y to X.
    
    Args:
        X: Reference embeddings (n_samples, n_dim)
        Y: Embeddings to align (n_samples, n_dim)
        return_disparity: Whether to return the final disparity
        max_iter: Maximum number of iterations
        tol: Convergence tolerance
    
    Returns:
        Aligned embeddings and optionally the final disparity
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    
    if X.shape != Y.shape:
        raise ValueError(f"Shape mismatch: X {X.shape} vs Y {Y.shape}")
    
    logger.info(f"🔄 Starting hyperbolic Procrustes alignment")
    logger.info(f"   X shape: {X.shape}, Y shape: {Y.shape}")
    
    # Initialize
    Y_aligned = Y.copy()
    prev_disparity = np.inf
    
    # Use origin as base point for tangent space
    base_point = np.zeros(X.shape[1])
    
    for iteration in range(max_iter):
        # Step 1: Log-map to tangent space
        X_tangent = log_map(X, base_point)
        Y_tangent = log_map(Y_aligned, base_point)
        
        # Step 2: Center the embeddings
        X_centered = center_embeddings(X_tangent)
        Y_centered = center_embeddings(Y_tangent)
        
        # Step 3: Apply Kabsch rotation
        R = kabsch_rotation(Y_centered, X_centered)
        Y_rotated = Y_centered @ R.T
        
        # Step 4: Exp-map back to hyperbolic space
        Y_aligned = exp_map(Y_rotated, base_point)
        
        # Check convergence
        if return_disparity:
            current_disparity = np.mean(hyperbolic_distance(X, Y_aligned))
            
            if abs(current_disparity - prev_disparity) < tol:
                logger.info(f"   ✅ Converged at iteration {iteration + 1}")
                logger.info(f"   Final disparity: {current_disparity:.6f}")
                break
            
            prev_disparity = current_disparity
            
            if iteration % 10 == 0:
                logger.info(f"   Iteration {iteration + 1}: disparity = {current_disparity:.6f}")
    
    if return_disparity:
        final_disparity = np.mean(hyperbolic_distance(X, Y_aligned))
        return Y_aligned, final_disparity
    else:
        return Y_aligned

def align_multiple_embeddings(embeddings_dict: dict, 
                            reference_key: str = None,
                            return_disparities: bool = False) -> dict:
    """
    Align multiple embeddings to a reference embedding.
    
    Args:
        embeddings_dict: Dictionary of embeddings {name: embeddings}
        reference_key: Key of reference embedding (if None, use first key)
        return_disparities: Whether to return disparities
    
    Returns:
        Dictionary of aligned embeddings and optionally disparities
    """
    if reference_key is None:
        reference_key = list(embeddings_dict.keys())[0]
    
    if reference_key not in embeddings_dict:
        raise ValueError(f"Reference key '{reference_key}' not found in embeddings_dict")
    
    reference = embeddings_dict[reference_key]
    aligned_embeddings = {reference_key: reference}
    disparities = {}
    
    for key, embeddings in embeddings_dict.items():
        if key == reference_key:
            continue
        
        logger.info(f"🔄 Aligning {key} to {reference_key}")
        
        if return_disparities:
            aligned_emb, disparity = hyperbolic_procrustes(
                reference, embeddings, return_disparity=True
            )
            aligned_embeddings[key] = aligned_emb
            disparities[key] = disparity
        else:
            aligned_emb = hyperbolic_procrustes(reference, embeddings)
            aligned_embeddings[key] = aligned_emb
    
    if return_disparities:
        return aligned_embeddings, disparities
    else:
        return aligned_embeddings

# Convenience function for backward compatibility
def procrustes_align(X: np.ndarray, Y: np.ndarray, 
                    return_disparity: bool = False) -> Tuple[np.ndarray, Optional[float]]:
    """
    Alias for hyperbolic_procrustes for backward compatibility.
    """
    return hyperbolic_procrustes(X, Y, return_disparity=return_disparity) 
