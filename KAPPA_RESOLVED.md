# KAPPA_RESOLVED.md

*Authoritative scientific record for the manuscript's κ claim.*
*Companion to KAPPA_COHERENCE.md. Written 2026-03-12 after completing Phases 1–3.*

---

## 1. The Dilemma (see KAPPA_COHERENCE.md for full narrative)

The manuscript claims κ = 1.247 ± 0.003 is *discovered* by neural networks and
agrees with Manning's theoretical prediction κ = (h ln 2)² ≈ 1.23 at h = 1.6 bits/nt.

KAPPA_COHERENCE.md established that this claim rests on three layers of circular or
fabricated evidence:

1. `five_seed_convergence.yaml` — fabricated before experiments ran; all κ values
   are `softplus_inverse(1.5) = 1.2475`, not from training.
2. `compact2` hardcodes κ = 1.25 at model creation via geoopt's `c_init` mechanism.
3. The aspirational YAMLs for viral and phylogenetic validation are future-dated
   theoretical predictions, not measurements.

---

## 2. The Telescope Experiment

**Method**: Joint LBFGS optimization of Poincaré coordinates + κ on GTDB patristic
distance matrices. Model-independent: no neural network involved. Uses the same
`HyperbolicTreeEmbedder` that recovered κ_nucleotide ≈ 1.23 for protein family trees
in the earlier protein_kappa experiment.

**Data**: GTDB r207 representative genomes, `build_patristic_matrix.py`.
- Bacteria (bac120 marker tree): N = 150 taxa
- Archaea (ar53 marker tree): N = 100 taxa
- Distance range: bacteria [0.017, 3.80] subs/site, archaea [0.032, 4.78] subs/site
- **Cross-domain distances are a constant sentinel (7.163) — combined matrix unusable.**
  Domain-split matrices required and used.

**Code**: `fit_kappa_phase1.py` (= `fit_kappa_telescope.py` from repo), run
2026-03-12 on inference machine (rohit@100.86.142.125).

---

## 3. The Result

| Domain   | N   | κ*     | 95% CI          | Stress  | 1.2475 ∈ CI? |
|----------|-----|--------|-----------------|---------|--------------|
| Bacteria | 150 | 17.63  | [17.43, 19.68]  | 0.0326  | **NO**       |
| Archaea  | 100 | 18.24  | [18.11, 19.97]  | 0.0505  | **NO**       |

- All three restarts (κ_init = 1.0, 5.0, 10.0) converged to κ* ≈ 17–18.
- Bootstraps are tight and consistent (std ≈ 0.6).
- Stress values are low — the embedding is genuinely good at this κ.
- κ = 1.2475 is **14× below** the measured optimal curvature.

---

## 4. The Entropy Chain (Phase 3)

**Measured from tokenized training data** (300 randomly sampled genomes):

| Quantity | Value |
|---|---|
| Mean nt per BPE token | **4.35 nt/token** |
| MLM loss at step 7000 (grouped_kappa_run) | **5.845 nats** |
| h_empirical = (mlm_loss / ln2) / mean_nt | **1.94 bits/nt** |
| κ_Manning_empirical = (h × ln2)² | **1.81** |
| Unigram ceiling (model-free) | h = 2.35 bits/nt → κ = 2.65 |

The entropy is h ≈ 1.94 bits/nt — **higher than the assumed 1.6 bits/nt** in the
manuscript. Manning's formula at the actual empirical h gives κ_Manning ≈ 1.81, not 1.25.

For κ = 1.2475 to be the Manning curvature, h would need to be exactly 1.611 bits/nt.
The model's MLM loss does not support this.

---

## 5. The Phase 2 Control (compact2 Telescope)

compact2 embeddings vs GTDB patristic distances (N = 250, mixed bacteria+archaea):

- Best Spearman ρ = **0.139** (at fixed κ = 1.25)
- Fitted κ wildly init-dependent: 0.20, 0.33, 0.48, 0.70, 2.31, 3.54
- No stable minimum; no convergence signal

**Interpretation**: compact2's embeddings at κ = 1.25 explain ~2% of variance
(ρ² ≈ 0.019) in GTDB patristic distances. The model is not encoding phylogenetic
structure at this curvature.

---

## 6. Neural Training Confirms the Artifact (Phase 4, grouped_kappa_run)

Training with learnable κ on 35,855 genomes, 900 genera, InfoNCE + MLM:

- Final κ = **0.775** at step 7000
- MLM = 5.845 nats, hex = 2.288 nats, dist_loss = 0.000 throughout

κ converged to ~0.78, not 1.247. The distance loss was identically zero throughout
training — the margin loss never fired, meaning the InfoNCE temperature was the
only signal on κ, and it pulled κ away from 1.25, not toward it.

This is the fourth independent experiment confirming KAPPA_COHERENCE.md's conclusion:
neural training does NOT converge to κ = 1.247.

---

## 7. The Unit Mismatch — Why the Telescope Gives κ* ≈ 18

The telescope measures the geometric curvature of the empirical tree **in the units of
the patristic distance metric** (GTDB marker gene substitutions per aligned site).

Manning's formula κ = (h ln 2)² was derived for **raw nucleotide positions** (bits/nt).

These are different metrics. The ratio between maximum and minimum patristic distances
in the GTDB bacterial tree is **224:1** after median normalization. A tree with that
depth-to-breadth ratio requires κ ≈ 18 to embed in 2D H² with normalized distances.
The protein family trees in the earlier experiment had ratios of ~10–50:1, which
is why they yielded κ* ≈ 1.2–1.3 — coincidentally near Manning's prediction, but
for reasons of tree size, not a universal law.

**The mapping between Manning's κ and tree-embedding κ requires:**

```
κ_tree = κ_Manning × (patristic_scale / nucleotide_scale)²
```

These scales differ by orders of magnitude at whole-tree depth. They would only agree
at a specific scale where patristic distances are calibrated in nucleotide information
units (bits/position), which GTDB distances are not.

---

## 8. The Verdict on κ = 1.2475

### What κ = 1.2475 is:
- `softplus_inverse(1.5)` — the geoopt default curvature initialization
- A round-number design choice (`c_init = 1.5`) that happens to sit within 3 significant
  figures of the Manning prediction at the assumed h = 1.6 bits/nt
- A design parameter that *works* for compact2 (coherent geometry, meaningful
  downstream tasks), analogous to fixing a coupling constant from theory

### What κ = 1.2475 is not:
- **Discovered by neural networks**: training converges to ~0.78
- **Matching the tree of life**: GTDB trees require κ* ≈ 17–18 in patristic units
- **Validated by five independent seeds**: those values are fabricated
- **Confirmed by viral experiments**: that YAML is aspirational, pre-execution
- **Confirmed by phylogenetic experiments**: those values are aspirational, pre-execution

### What the experiments support:
- κ = 1.25 works as a *fixed design parameter* for a genome embedding model
- The geometry is coherent: E1–E2 validation (cos sim ≈ 0.977, cross-domain transfer ≈ 0.96)
- The claim that this specific value is *discovered* or *forced by information theory*
  is not supported by any of the measurements

---

## 9. What compact2's Success Means

compact2 with hardcoded κ = 1.25 produces a coherent latent space that passes E1–E2.
This is genuine and reproducible. But it does not mean 1.25 is the uniquely correct
curvature. It means:

- A Poincaré ball at κ = 1.25 is a reasonable geometry for organizing 47k genomes
- The model learned useful representations within this geometry
- The specific value 1.25 was chosen by the researcher (via `c_init = 1.5`), not
  discovered by the model
- The downstream tasks (cross-domain transfer, E2 structure) are real, regardless of
  whether 1.25 is special

The analogy to "fixing coupling constants from theory and validating the detector"
holds only if the theory genuinely predicted 1.25. The measurement shows it did not.

---

## 10. Path Forward

The manuscript requires revision in the following respects:

1. **Remove the claim that neural networks *discover* κ = 1.247**. Replace with:
   "We train with κ = 1.25 as a fixed design choice motivated by Manning's formula at
   the assumed h = 1.6 bits/nt. Validation of whether this value is empirically
   optimal requires a telescope experiment with distances in nucleotide information units."

2. **Remove the five_seed_convergence.yaml evidence.** Mark that file as fabricated
   (done in Phase 0).

3. **Add the Phase 1 telescope result**: κ*_bacteria = 17.6, κ*_archaea = 18.2 from
   GTDB patristic distances. State clearly that these are in substitution/site units
   and are not directly comparable to Manning's κ.

4. **Add the Phase 3 entropy result**: h_empirical ≈ 1.94 bits/nt from the model's
   MLM loss; κ_Manning_empirical ≈ 1.81. The assumed h = 1.6 bits/nt appears too low.

5. **The canonical telescope for the manuscript's claim** would require patristic
   distances expressed in raw nucleotide mutual information (bits), not marker gene
   substitutions. This experiment has not yet been run and requires a different
   distance computation.

---

## Appendix: Summary of All κ Measurements

| Experiment | κ measured | Method | Date | Status |
|---|---|---|---|---|
| five_seed_convergence.yaml | 1.247 | "5-seed training" | 2025-01-01 | **FABRICATED** |
| compact2 model | 1.250 | `c_init=1.5` at model creation | pre-2026 | Design choice |
| grouped_kappa_run step 7000 | 0.775 | Learnable κ, InfoNCE | 2026-03-12 | Real |
| Phase 2 telescope (compact2) | init-dependent 0.2–3.5 | κ fit to GTDB patristic | 2026-03-12 | Real, no signal |
| Phase 1 telescope (bacteria) | **17.63** | Pure tree LBFGS | 2026-03-12 | Real |
| Phase 1 telescope (archaea) | **18.24** | Pure tree LBFGS | 2026-03-12 | Real |
| Phase 3: κ_Manning from MLM | **1.81** | h_empirical from loss | 2026-03-12 | Real |
| Manning theory (h=1.6 bits/nt) | 1.228 | (h ln 2)² | — | Theoretical |
| Manning theory (h=1.94 bits/nt) | 1.813 | (h ln 2)² | — | Consistent w/ Phase 3 |
