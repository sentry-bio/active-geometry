# KAPPA_PHASE2_DIAGNOSTIC.md

*Phase 2 angular refinement comparison: κ=1.25 vs κ=1.39*
*Written 2026-03-13 after stopping both runs at 21–24 epochs.*
*Raw logs: `phase2_kappa125_train.log`, `phase2_kappa139_train.log`*

---

## 1. Run Configuration

| Parameter | κ=1.25 | κ=1.39 |
|-----------|--------|--------|
| Script | `train_v15_5_phase2_kappa125.py` | `train_v15_5_phase2_kappa139.py` |
| Checkpoint source | `phase1_kappa125/best.pt` | `retrain_kappa139_scratch/best.pt` |
| Epochs run | 21 (stopped) | 24 (completed) |
| Radial head | FROZEN (a=1.191, b=0.087) | FROZEN (a=1.272, b=0.156) |
| kappa_prior_weight | 0.0 (no prior) | 0.0 (no prior) |
| AngSep weight/threshold | 0.15 / 0.85 | 0.15 / 0.85 |
| Batch size / grad accum | 10 / 14 | 10 / 14 |

**Radial state entering Phase 2** is itself diagnostic: κ=1.39 has b=0.156 (78% higher than κ=1.25's b=0.087). A larger b correction implies the radial path was doing more compensatory work during Phase 1 — the geometry at 1.39 required more post-hoc radius adjustment to maintain monotonic ordering.

---

## 2. κ Drift Table (epoch-by-epoch)

Both runs nominally "froze" κ, yet both drifted upward throughout. The drift is a known artifact: `requires_grad_(False)` is applied after the optimizer is constructed, so Adam momentum continues carrying κ forward.

### κ=1.25 run

| Epoch | Sampler | κ | Δκ from init |
|-------|---------|------|--------------|
| init (checkpoint load) | — | 1.2506 | — |
| 1 | UF | 1.2521 | +0.0015 |
| 3 | FS | 1.2543 | +0.0037 |
| 5 | UF | 1.2589 | +0.0083 |
| 6 | FS | 1.2591 | +0.0085 |
| 9 | FS | 1.2639 | +0.0133 |
| 12 | FS | 1.2687 | +0.0181 |
| 15 | FS | 1.2738 | +0.0232 |
| 18 | FS | 1.2787 | +0.0281 |
| 21 | FS | 1.2837 | +0.0331 |

**Drift rate: +0.00158/epoch** — slow, nearly uniform.

### κ=1.39 run

| Epoch | Sampler | κ | Δκ from script reset (1.39) |
|-------|---------|------|------------------------------|
| init (checkpoint load → reset to 1.39) | — | 1.3900 | — |
| 1 | UF | 1.3921 | +0.0021 |
| 3 | FS | 1.3945 | +0.0045 |
| 5 | UF | 1.3991 | +0.0091 |
| 6 | FS | 1.3993 | +0.0093 |
| 9 | FS | 1.4041 | +0.0141 |
| 12 | FS | 1.4091 | +0.0191 |
| 15 | FS | 1.4141 | +0.0241 |
| 18 | FS | 1.4191 | +0.0291 |
| 21 | FS | 1.4242 | +0.0342 |
| 24 | FS | 1.4291 | +0.0391 |

**Drift rate: +0.00154/epoch** — almost identical to κ=1.25.

### Interpretation

The drift rates are **indistinguishable** (0.00158 vs 0.00154/epoch). This rules out any interpretation of the drift as the geometry seeking an attractor. It is pure optimizer momentum artifact — the Adam state was initialized during Phase 1 with κ trainable, and the accumulated first/second moments keep pushing κ upward at a constant rate regardless of where we "freeze" it. Neither 1.25 nor 1.39 shows the drift decelerating, which would be the signature of approaching a true fixed point.

**Crucially**: if either value were the thermodynamic fixed point, we would expect the gradient from the geometry to oppose the drift. Neither does. The drift is an artifact, not signal.

---

## 3. FS Comb Trajectory (full-shuffle epochs only)

These are the only fair performance comparisons — UF epochs see different data distributions.

| FS Epoch | κ=1.25 Comb | κ=1.39 Comb | Δ (1.25 − 1.39) |
|----------|-------------|-------------|-----------------|
| Baseline (entering Phase 2) | 19.7 | 23.2 | −3.5 |
| Ep 3 | 26.3 | 26.8 | −0.5 |
| Ep 6 | 25.9 | 27.3 | −1.4 |
| Ep 9 | 27.1 | 26.9 | +0.2 |
| Ep 12 | 27.0 | 26.7 | +0.3 |
| Ep 15 | **28.9** | 28.6 | +0.3 |
| Ep 18 | 28.7 | **28.8** | −0.1 |
| Ep 21 | 28.2 | 28.1 | +0.1 |
| Ep 24 | — | 27.6 | — |

**Peak:** κ=1.25 → 28.93 at ep15. κ=1.39 → 28.80 at ep18. Difference: **0.13 Comb** — within noise.

**Phase 2 improvement** (peak − baseline entry):
- κ=1.25: 28.93 − 19.7 = **+9.2**
- κ=1.39: 28.80 − 23.2 = **+5.6**

κ=1.25 gained significantly more through Phase 2, closing the gap opened during Phase 1. By ep9, the two runs are statistically indistinguishable despite κ=1.25 entering Phase 2 with a 3.5 Comb deficit.

Both plateau in the **28–29 Comb band** and show no sign of escaping it through continued training. This plateau is discussed in Section 6.

---

## 4. UF/FS Oscillation Pattern

UF (UniformFamily) samples one example per family from 1664 groups — the hardest retrieval setting. FS (full-shuffle) samples uniformly from 37k genomes.

### κ=1.25 UF Comb

`1.7 → 2.1 → 1.9 → 2.0 → 2.1 → 2.2 → 2.2 → 2.1 → 2.5 → 2.4 → 2.4 → 2.7 → 2.3 → 2.8`

### κ=1.39 UF Comb

`2.1 → 2.4 → 1.8 → 2.1 → 2.3 → 2.3 → 2.2 → 2.2 → 2.3 → 2.7 → 2.4 → 3.2 → 2.5 → 3.2`

Both are locked in the **2–3 range** for the entire Phase 2. Neither trajectory shows any resolution of the UF/FS gap. κ=1.39 reaches slightly higher UF peaks (3.2 at ep17/20 vs 2.8 at ep20 for κ=1.25) but neither is trending upward.

The predicted "smoking gun" — UF Comb beginning to track FS Comb as prototypes crystallize — did not materialise in either run. This is addressed in Section 6.

---

## 5. GeoProbe Deep Comparison

Full GeoProbe readings at the baseline entry and at selected FS epochs.

### 5a. Baseline (entering Phase 2)

| Metric | κ=1.25 | κ=1.39 |
|--------|--------|--------|
| Arch retrieval | **30.0%** | 23.4% |
| hard_mono | 0.382 | 0.254 |
| δ (quartet separation) | 0.076 | 0.082 |
| proto_entropy B/A/E | 0.715/0.738/0.681 | 0.775/0.708/0.774 |
| angular_spread B/A/E | **0.770/0.762/0.746** | 0.692/0.695/0.702 |

κ=1.25 enters Phase 2 with meaningfully better arch retrieval (30% vs 23.4%) and angular spread. κ=1.39 has slightly higher proto_entropy overall (prototypes less specialised, more uniform). The angular_spread advantage for κ=1.25 is notable — the ball at κ=1.25 (radius ≈ 0.894) is slightly wider, giving prototypes more room to spread angularly before being compressed toward the boundary.

### 5b. AngSep activation

| Run | Ep1 | Ep6 | Ep12 | Ep18 | Ep21 |
|-----|-----|-----|------|------|------|
| κ=1.25 | 0.0004 | 0.0017 | 0.0024 | 0.0026 | 0.0029 |
| κ=1.39 | 0.0008 | 0.0021 | 0.0028 | 0.0031 | 0.0033 |

**Neither run meaningfully activated the angular separation loss.** The threshold is 0.85; angular_spread across both runs sits at 0.69–0.77. The AngSep loss fires only for the very small fraction of prototype pairs that are dangerously close (cos sim > 0.85). This is a structural ceiling: Phase 2 with threshold=0.85 cannot drive meaningful prototype spreading when the geometry is already near that spread level. Phase 2b (threshold lowering) was always the intended next step.

### 5c. hard_mono trajectory (FS epochs only)

| FS Epoch | κ=1.25 | κ=1.39 |
|----------|--------|--------|
| Baseline | 0.382 | 0.254 |
| Ep 3 | 0.318 | 0.512 |
| Ep 6 | 0.495 | 0.335 |
| Ep 9 | 0.324 | 0.524 |
| Ep 12 | 0.492 | 0.498 |
| Ep 15 | 0.328 | **0.689** |
| Ep 18 | 0.341 | **0.691** |
| Ep 21 | **0.509** | 0.501 |

This is the most striking divergence. κ=1.39 achieves hard_mono peaks of **0.689–0.691** at ep15–18, while κ=1.25 peaks at 0.509. However, three observations qualify this:

1. Both show the same UF/FS alternation pattern — hard_mono is high on FS epochs and drops on UF. The "peak" is a FS measurement, not a stable geometric property.
2. κ=1.39's spike at ep15–18 is followed by a drop back to 0.501 at ep21/24, suggesting it is not stabilising at that level.
3. κ=1.25's hard_mono enters Phase 2 at 0.382 (already higher than κ=1.39's 0.254 baseline), reflecting better Phase 1 radial organisation.

The hard_mono advantage for κ=1.39 at ep15–18 is real but transient. It reflects a brief period where the radial ordering (from the higher b=0.156) aligns well with the FS sampling pattern.

### 5d. Arch retrieval trajectory

| FS Epoch | κ=1.25 | κ=1.39 |
|----------|--------|--------|
| Baseline | **30.0%** | 23.4% |
| Ep 3 | 38.8% | **41.6%** |
| Ep 6 | **40.6%** | 38.2% |
| Ep 9 | 33.8% | 23.6% |
| Ep 12 | 37.6% | 34.2% |
| Ep 15 | 33.2% | 26.2% |
| Ep 18 | 34.4% | **39.2%** |
| Ep 21 | 37.2% | 33.8% |
| Ep 24 | — | 34.6% |

Both oscillate substantially (15–20% swings). Neither shows a rising trend. κ=1.25 has the higher baseline (30% vs 23.4%) — a genuine Phase 1 advantage. Peak values are comparable (40.6% vs 41.6%).

### 5e. 2D_var (GCAS — fraction of embedding variance explained by first 2 dims)

| GCAS Epoch | κ=1.25 | κ=1.39 |
|------------|--------|--------|
| Ep 5 | 13.7% | 14.8% |
| Ep 10 | 15.0% | 16.0% |
| Ep 15 | 17.2% | **17.9%** |
| Ep 20 | 16.5% | 17.5% |

κ=1.39 consistently runs ~1% higher in 2D_var. This is consistent with a slightly more compact angular arrangement — at higher curvature the ball is smaller, so 128D embeddings project onto fewer effective dimensions. This is not a clear advantage; it may simply reflect that κ=1.39's geometry is more compressed.

---

## 6. The Phase 2 Ceiling: Why Both Plateau at ~28-29

Both runs plateau in the same 28–29 Comb band and fail to resolve the UF/FS oscillation. This was not predicted to be symmetric — the earlier hypothesis was that κ=1.25 would break through 28-29 while κ=1.39 would not. That hypothesis was wrong. The symmetry itself is informative.

**The bottleneck is not curvature.** It is the AngSep activation threshold. Phase 2 with threshold=0.85 and weight=0.15 is structurally incapable of driving the angular separation needed to lift UF Comb. The AngSep gradient is ≈0 for 95%+ of prototype pairs at either curvature. Without angular separation pressure, the ODE flow can rearrange angular positions locally (improving FS Comb through better family clustering) but cannot force the global disambiguation needed for UF retrieval.

**What Phase 2b does**: lowers the AngSep threshold from 0.85 to ~0.70 and raises the weight, forcing separation across a much larger fraction of prototype pairs. This is precisely the missing mechanism. Both κ=1.25 and κ=1.39 would likely benefit similarly from Phase 2b — the ceiling is architectural, not geometric.

**What this means for κ comparison**: the curvature competition cannot be settled at this phase of training. Phase 2 shows they are roughly equivalent. The tiebreakers are:
- κ=1.25 started Phase 2 with higher arch retrieval and angular spread
- κ=1.25 gained more through Phase 2 (+9.2 vs +5.6 Comb)
- κ=1.39 showed better transient hard_mono peaks
- Both are at the same ceiling now

---

## 7. Summary Scoreboard

| Metric | κ=1.25 | κ=1.39 | Edge |
|--------|--------|--------|------|
| Phase 1 best Comb | 21.9 | 23.1 | 1.39 |
| Phase 2 baseline Comb | 19.7 | 23.2 | **1.39** |
| Phase 2 peak Comb | **28.93** | 28.80 | 1.25 (marginal) |
| Phase 2 Comb gain | **+9.2** | +5.6 | **1.25** |
| Baseline arch retrieval | **30.0%** | 23.4% | **1.25** |
| Peak arch retrieval | 40.6% | 41.6% | tie |
| Angular spread (baseline) | **0.770/0.762/0.746** | 0.692/0.695/0.702 | **1.25** |
| Peak hard_mono | 0.509 | **0.689** | 1.39 |
| UF Comb peak | 2.8 | **3.2** | 1.39 |
| UF/FS oscillation resolved | No | No | tie |
| AngSep activation | 0.0029 | 0.0034 | tie (both ≈0) |
| κ drift rate | 0.00158/ep | 0.00154/ep | tie (artifact) |
| Radial b entering P2 | 0.087 | 0.156 | — |

**Overall**: indistinguishable at this stage. The Phase 2 ceiling is shared and is not a function of κ.

---

## 8. What the Geometry Looks Like Right Now

At ep21 for κ=1.25:
- `r mean=1.74, std=0.090, range=[0.99, 1.81]` — embeddings tightly clustered in radius, not yet exploiting the ball depth
- `corr(r,r0)=0.999` — perfect rank-order preservation of Phase 1 radial targets; Phase 2 has not disturbed radial structure
- `proto_std B:0.233 A:0.304 E:0.416` — eukaryote prototypes most spread (443 families competing for angular positions), archaea most concentrated
- `proto_entropy B:0.384 A:0.341 E:0.574` — significant collapse from baseline (0.715/0.738/0.681); prototypes have specialised toward specific families, especially archaea

The geometry is **crystallising but not crystallised**. Prototypes have found family-specific angular niches (collapsing entropy), the radial ordering is intact, and arch retrieval is sitting at 30–40%. The Phase 2c annealing step — which unfreezes the radial head and trains the full model simultaneously — is the mechanism that historically (compact2: Comb→64.6) pushes this half-crystallised state into a fully organised manifold.

---

## 9. Recommended Next Step

Phase 2c (compact2 style) for **both** checkpoints:
- Unfreeze radial head (b can grow beyond its Phase 1 value)
- Train encoder + ODE + heads + radial simultaneously
- Full 50k manifest
- ~30 epochs

The Phase 2c run on κ=1.25's current best.pt is the highest-leverage experiment. If κ=1.25's Comb climbs toward compact2's 64.6, that validates the full training pipeline and gives a new baseline checkpoint. If it plateaus well below 64.6, that points to a Phase 2b step (AngSep threshold lowering) being needed before 2c.

For a principled head-to-head, run Phase 2c on both κ=1.25 (ep21 best.pt) and κ=1.39 (ep24 best.pt) with identical hyperparameters and compare final Comb, E1/E2 validation scores, and cross-domain transfer matrices.
