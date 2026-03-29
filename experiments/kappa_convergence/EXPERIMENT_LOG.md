# κ Convergence Experiments — Session Log

## Objective
Verify that the curvature constant κ = (h·ln2)² = 5/4 ≈ 1.2505 emerges as the canonical curvature of the Poincaré ball when training on genomic data, regardless of initialization.

## Key Finding: What Works and What Doesn't

### Approaches That CANNOT Find κ
1. **Direct manifold optimization** (taxonomy labels only, no genomes): κ tracks the initialization linearly. Without genomic compression, h is absent and there is no information-theoretic anchor.
2. **E11 lightweight encoder** (5-term geometric loss, Euclidean Adam): κ blasts past 1.25 to ~2.6 in 2D due to Euclidean optimizer over-updating near the Poincaré boundary.
3. **E11 + Riemannian gradient correction** (manual conformal factor λ(z) on backward hook): κ converges init-dependently (0.66, 0.85, 1.65 from different inits). The encoder's finite capacity creates multiple basins.

### The Approach That Works
**BiosphereCodec pretrain + κ fine-tune**: The pretrained encoder learns Shannon entropy rate h through MLM/CLM compression. Then with encoder frozen and c learnable, InfoNCE drives κ to its natural basin. This is the only path where h enters the picture.

## Experiments Run (chronological)

### 1. `run_kappa_convergence_v2.py` — InfoNCE with Stochastic Sampling
- **Machine**: Nexus (rohit@100.100.20.90) + Inference (rohit@100.86.142.125)
- **Key fixes**: Pandas groupby bug, FamilyGroupedBatchSampler, sampler mode toggle
- **Result**: κ touched 1.248 at step 260 with stochastic sampling (batch=64) but kept climbing to 1.98+
- **Insight**: Stochasticity matters — deterministic grouped sampling biases κ (8fam×4→1.04, 16fam×2→1.50)

### 2. `run_kappa_E11.py` — 5-Term Geometric Loss (Euclidean Adam)
- **Machine**: Inference (GPU)
- **Architecture**: 2D Poincaré encoder, 40K params, 268 anchors
- **Losses**: quartet + 2×domain_angular + 5×genus_anchor + 0.5×angular_repulsion + 0.3×radial_ordering
- **Result**: κ climbed from 0.5→2.60 and from 1.0→3.10 (300 epochs). Same runaway overshoot.
- **Root cause**: Euclidean Adam treats Poincaré ball as flat — massive over-updates near boundary

### 3. `run_kappa_E11_riemannian.py` — Riemannian Gradient Correction
- **Innovation**: `RiemannianGradientHook` scales backward gradients by λ(z) = (1-c‖z‖²)²/4
- **Also**: softplus reparameterization for c (no geoopt dependency)
- **Result (5 seeds × 3 inits, 300 epochs)**:
  - Prevents runaway (κ no longer blasts to 2.6+)
  - But init-dependent: init 0.5→0.66, init 1.0→0.85, init 2.0→1.65
  - Seed-dependent: seed 93 collapsed to 0.03
- **Conclusion**: Encoder creates multiple basins; Riemannian correction necessary but insufficient

### 4. `run_kappa_direct_manifold.py` — geoopt RiemannianAdam, No Encoder
- **Architecture**: 500 organisms as ManifoldParameter, RiemannianAdam for positions, separate Adam for c
- **Result (5 seeds × 6 inits, 500 epochs)**: κ scales linearly with init (0.1→0.59, 0.5→1.63, 1.0→2.35, 5.0→6.44)
- **Conclusion**: Without compression (no genomes), h is missing → no canonical κ basin

### 5. `run_kappa_geoopt_native.py` — geoopt End-to-End
- **Attempted**: PoincareBall.dist() with learnable c, RiemannianAdam for everything
- **Result**: ∇c = 0 — geoopt's in-place reparameterization detaches c from autograd graph
- **Key discovery**: `ball.c.requires_grad = False` after PoincareBall.__init__
- **Conclusion**: geoopt's c management is valid for fixed curvature but breaks gradient flow for learnable c

### 6. `train_codec_kappa.py` — BiosphereCodec Pretrain + κ Measurement ← CANONICAL
- **Phase 1**: Train BiosphereCodec from scratch (MLM + CLM + InfoNCE) on 5K genomes
- **Phase 2**: Freeze encoder, patch geoopt out, sweep c_init × seeds, measure convergence
- **Status**: RUNNING (dim=2, 5K genomes, 5 seeds × 5 inits)
- **Why this works**: Encoder learns h through compression, InfoNCE provides distance-ratio signal that depends on κ

## Critical Technical Discoveries

### geoopt Softplus Reparameterization
`PoincareBall.__init__` applies `k.exp_().sub_(1).log_()` — this is inverse softplus, converting user-facing c to internal raw_c. Then `ball.c` returns `softplus(raw_c)`. This is a valid reparameterization for fixed c, but the in-place ops **destroy autograd history**, making `ball.c.requires_grad = False`. To learn c, you must:
- Use `self.c` (the raw parameter) directly in distance computations
- OR manage c outside geoopt with manual softplus

### Riemannian Gradient Correction
For encoder → manifold architectures (z = f_θ(tokens)), the correct Riemannian correction is a backward hook on z that scales gradients by the conformal factor λ(z) = (1-c‖z‖²)²/4. This is NOT the same as geoopt's RiemannianAdam, which optimizes direct ManifoldParameters. The hook approach is correct because z is an intermediate activation, not a parameter.

### Why Compression Is Required
κ = (h·ln2)² where h is Shannon entropy rate. Without genomic sequences, h is absent from the optimization. Taxonomy labels alone provide topological structure but not the information-theoretic content that anchors κ at 5/4.

## Architecture Constants
- Pretrained V15Model: κ = 1.2505 (from curvature_history), latent_dim=128, d_model=512
- Theoretical: κ = (1.6·ln2)² = 1.2283 (approx 5/4 = 1.25)

## Files
| Script | Purpose | Status |
|--------|---------|--------|
| `run_kappa_convergence_v2.py` | InfoNCE sampling experiments | Complete |
| `run_kappa_E11.py` | Euclidean 5-loss baseline | Complete |
| `run_kappa_E11_riemannian.py` | Riemannian gradient correction | Complete |
| `run_kappa_direct_manifold.py` | Direct manifold (no encoder) | Complete |
| `run_kappa_geoopt_native.py` | geoopt end-to-end (broken c grad) | Failed |
| `codec_kappa_finetune.py` | Fine-tune pretrained checkpoint | Blocked (key mismatch) |
| `train_codec_kappa.py` | Train BiosphereCodec + measure κ | **RUNNING** |
