/-
  The Metric Packing Theorem
  ==========================

  This file proves the counting core of the addressability limit.

  A retained representation supplies, at every generative depth R:

  * a finite, nonempty set of represented histories;
  * fixed-resolution separation between distinct represented histories;
  * containment in a ball of radius r(R).

  Mathlib's `Metric.packingNumber` is the maximal cardinality, in extended
  naturals, of a separated subset of a metric set. Hence the represented count
  is bounded by the packing number of its containing ball at every depth.
  Conversely, when that number is finite, an exact finite block code attaining
  it exists. Thus operational finite-block address capacity and metric packing
  capacity are identical in every proper metric host.

  If the following three finite rates exist:

      log represented-count / R  → β,
      r(R) / R                    → c,
      log packing-count / ρ       → h_pack,

  then

      β ≤ c · h_pack.

  The conclusion is exactly `Addressable β c h_pack` from
  `Addressability.lean`. No curvature, isotropy, saturation, dimension, or
  physical dynamics enters this proof. The block-achievability theorem does not
  claim that optimal codebooks across depths are nested, causal, or preserve a
  source hierarchy's relational metric.
-/

import ActiveGeometry.Addressability
import Mathlib.Topology.MetricSpace.CoveringNumbers
import Mathlib.Topology.MetricSpace.ProperSpace
import Mathlib.Topology.Order.OrderClosed

namespace ActiveGeometry.Packing

open Filter
open scoped ENNReal NNReal Topology

variable {M : Type*} [MetricSpace M]

/-- Packing numbers are finite at the chosen resolution on every
    non-negative-radius ball. This is the local finiteness needed to convert
    Mathlib's extended-natural packing number into an ordinary count. -/
def HasFinitePacking (o : M) (ε : ℝ≥0) : Prop :=
  ∀ ⦃ρ : ℝ⦄, 0 ≤ ρ →
    Metric.packingNumber ε (Metric.closedBall o ρ) ≠ ⊤

/-- The ordinary-natural packing count of a metric ball. When the extended
    packing number is infinite, `ENat.toNat` returns zero; all theorems using
    this definition therefore require `HasFinitePacking`. -/
noncomputable def packingCount (o : M) (ε : ℝ≥0) (ρ : ℝ) : ℕ :=
  (Metric.packingNumber ε (Metric.closedBall o ρ)).toNat

/-- Logarithmic packing growth per unit radius. -/
noncomputable def packingRate (o : M) (ε : ℝ≥0) (ρ : ℝ) : ℝ :=
  Real.log (packingCount o ε ρ : ℝ) / ρ

/-- A finite represented hierarchy at fixed metric resolution. Histories are
    represented by the points themselves; separation is the faithfulness
    condition, and containment records the addressing radius. -/
structure RetainedRepresentation (o : M) (ε : ℝ≥0) where
  points : ℕ → Finset M
  radius : ℕ → ℝ
  resolution_pos : 0 < ε
  points_nonempty : ∀ R, (points R).Nonempty
  /-- Retention: histories present at depth `R` persist to depth `R + 1`. This
      is the load-bearing meaning of the word "retained"; the bound itself does
      not need it, but it is what makes `β` a growth rate of retained history. -/
  points_monotone : ∀ R, points R ⊆ points (R + 1)
  separated : ∀ R, Metric.IsSeparated ε (points R : Set M)
  contained : ∀ R, (points R : Set M) ⊆ Metric.closedBall o (radius R)
  radius_pos : ∀ R, 0 < radius R
  radius_tendsto : Tendsto radius atTop atTop

/-- Retained-history growth observed at depth `R`. -/
noncomputable def representedRate
    {o : M} {ε : ℝ≥0} (rep : RetainedRepresentation o ε) (R : ℕ) : ℝ :=
  Real.log ((rep.points R).card : ℝ) / (R : ℝ)

/-- Radial distance used per generative step at depth `R`. -/
noncomputable def radialRate
    {o : M} {ε : ℝ≥0} (rep : RetainedRepresentation o ε) (R : ℕ) : ℝ :=
  rep.radius R / (R : ℝ)

/-- Every finite separated set contained in a ball has cardinality at most the
    exact packing number of that ball. -/
theorem card_le_packingCount
    (o : M) (ε : ℝ≥0) (ρ : ℝ) (s : Finset M)
    (hρ : 0 ≤ ρ) (hfinite : HasFinitePacking o ε)
    (hsep : Metric.IsSeparated ε (s : Set M))
    (hcontained : (s : Set M) ⊆ Metric.closedBall o ρ) :
    s.card ≤ packingCount o ε ρ := by
  have henat :
      ((s.card : ℕ) : ℕ∞) ≤
        Metric.packingNumber ε (Metric.closedBall o ρ) := by
    simpa using hsep.encard_le_packingNumber hcontained
  have htop :
      Metric.packingNumber ε (Metric.closedBall o ρ) ≠ ⊤ :=
    hfinite hρ
  simpa [packingCount] using ENat.toNat_le_toNat henat htop

/-- Finite-block achievability: whenever the ball has finite packing number,
    there is an actual finite codebook whose cardinality equals it. Together
    with `card_le_packingCount`, this identifies operational block-address
    capacity exactly with metric packing capacity. No tree, hyperbolicity, or
    saturation hypothesis is involved.

    This is deliberately a block result: codebooks at different radii need not
    be nested, causal, or preserve a source's relational metric. -/
theorem exists_optimal_blockCode
    (o : M) (ε : ℝ≥0) (ρ : ℝ)
    (hρ : 0 ≤ ρ) (hfinite : HasFinitePacking o ε) :
    ∃ s : Finset M,
      Metric.IsSeparated ε (s : Set M) ∧
      (s : Set M) ⊆ Metric.closedBall o ρ ∧
      s.card = packingCount o ε ρ := by
  have htop :
      Metric.packingNumber ε (Metric.closedBall o ρ) ≠ ⊤ :=
    hfinite hρ
  obtain ⟨C, hcontained, hCfinite, hsep, hcard⟩ :=
    Metric.exists_set_encard_eq_packingNumber htop
  refine ⟨hCfinite.toFinset, ?_, ?_, ?_⟩
  · simpa only [hCfinite.coe_toFinset] using hsep
  · simpa only [hCfinite.coe_toFinset] using hcontained
  · rw [← Set.ncard_eq_toFinset_card C hCfinite]
    simpa only [Set.ncard_def, packingCount] using congrArg ENat.toNat hcard

/-- Finite-depth addressability: represented histories cannot outnumber the
    exact packing count of their containing ball. -/
theorem represented_card_le_packingCount
    (o : M) (ε : ℝ≥0) (rep : RetainedRepresentation o ε)
    (hfinite : HasFinitePacking o ε) (R : ℕ) :
    (rep.points R).card ≤ packingCount o ε (rep.radius R) :=
  card_le_packingCount o ε (rep.radius R) (rep.points R)
    (rep.radius_pos R).le hfinite (rep.separated R) (rep.contained R)

/-- The finite-depth count inequality becomes a pointwise inequality between
    normalized logarithmic rates away from the irrelevant depth `R = 0`. -/
theorem representedRate_le_capacity_eventually
    (o : M) (ε : ℝ≥0) (rep : RetainedRepresentation o ε)
    (hfinite : HasFinitePacking o ε) :
    ∀ᶠ R in atTop,
      representedRate rep R ≤
        radialRate rep R * packingRate o ε (rep.radius R) := by
  filter_upwards [eventually_gt_atTop 0] with R hR
  have hcount := represented_card_le_packingCount o ε rep hfinite R
  have hcard_pos : 0 < ((rep.points R).card : ℝ) := by
    exact_mod_cast Finset.card_pos.mpr (rep.points_nonempty R)
  have hcount_real :
      ((rep.points R).card : ℝ) ≤
        (packingCount o ε (rep.radius R) : ℝ) := by
    exact_mod_cast hcount
  have hlog :
      Real.log ((rep.points R).card : ℝ) ≤
        Real.log (packingCount o ε (rep.radius R) : ℝ) :=
    Real.log_le_log hcard_pos hcount_real
  have hR_real : 0 < (R : ℝ) := by exact_mod_cast hR
  have hr := rep.radius_pos R
  unfold representedRate radialRate packingRate
  calc
    Real.log ((rep.points R).card : ℝ) / (R : ℝ)
        ≤ Real.log (packingCount o ε (rep.radius R) : ℝ) / (R : ℝ) :=
      (div_le_div_iff_of_pos_right hR_real).2 hlog
    _ = (rep.radius R / (R : ℝ)) *
          (Real.log (packingCount o ε (rep.radius R) : ℝ) / rep.radius R) := by
      field_simp [hR_real.ne', hr.ne']

/-- The metric packing theorem. The three rate hypotheses are independent:
    retained-history growth, radial calibration, and host packing growth are
    not defined from one another. -/
theorem faithful_representation_addressable
    (o : M) (ε : ℝ≥0) (rep : RetainedRepresentation o ε)
    (hfinite : HasFinitePacking o ε)
    (β c hpack : ℝ)
    (hdemand : Tendsto (representedRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 hpack)) :
    Addressability.Addressable β c hpack := by
  have hsampled :
      Tendsto (fun R ↦ packingRate o ε (rep.radius R)) atTop (𝓝 hpack) := by
    simpa [Function.comp_def] using hpacking.comp rep.radius_tendsto
  apply le_of_tendsto_of_tendsto hdemand (hradial.mul hsampled)
  exact representedRate_le_capacity_eventually o ε rep hfinite

/-- Retention makes represented counts nondecreasing in depth. -/
theorem represented_card_mono {o : M} {ε : ℝ≥0}
    (rep : RetainedRepresentation o ε) (R : ℕ) :
    (rep.points R).card ≤ (rep.points (R + 1)).card :=
  Finset.card_le_card (rep.points_monotone R)

/-- Local finiteness of packing is a theorem, not an assumption, for the entire
    intended host class: proper metric spaces (ℝⁿ, hyperbolic space, and every
    complete Riemannian manifold via Hopf–Rinow). Closed balls are compact,
    hence totally bounded, hence finitely covered at resolution `ε / 2`, which
    bounds the packing number at resolution `ε`. This is where a positive
    resolution is genuinely required. -/
theorem hasFinitePacking_of_properSpace [ProperSpace M] (o : M) {ε : ℝ≥0}
    (hε : ε ≠ 0) : HasFinitePacking o ε := by
  intro ρ _hρ
  have hδ0 : ε / 2 ≠ 0 := div_ne_zero hε two_ne_zero
  have h2δ : 2 * (ε / 2) = ε := by
    rw [mul_comm]; exact div_mul_cancel₀ ε two_ne_zero
  have htb : TotallyBounded (Metric.closedBall o ρ) :=
    (isCompact_closedBall o ρ).totallyBounded
  obtain ⟨N, _hNsub, hNfin, hNcov⟩ :=
    Metric.exists_finite_isCover_of_totallyBounded hδ0 htb
  have hpack_le :
      Metric.packingNumber ε (Metric.closedBall o ρ)
        ≤ Metric.externalCoveringNumber (ε / 2) (Metric.closedBall o ρ) := by
    calc
      Metric.packingNumber ε (Metric.closedBall o ρ)
          = Metric.packingNumber (2 * (ε / 2)) (Metric.closedBall o ρ) := by
            rw [h2δ]
      _ ≤ Metric.externalCoveringNumber (ε / 2) (Metric.closedBall o ρ) :=
            Metric.packingNumber_two_mul_le_externalCoveringNumber (ε / 2) _
  have hcov_le :
      Metric.externalCoveringNumber (ε / 2) (Metric.closedBall o ρ) ≤ N.encard :=
    hNcov.externalCoveringNumber_le_encard
  exact (lt_of_le_of_lt (hpack_le.trans hcov_le) hNfin.encard_lt_top).ne

/-- The finite-block capacity identity with local finiteness discharged for a
    proper metric host. -/
theorem exists_optimal_blockCode_of_properSpace
    [ProperSpace M] (o : M) {ε : ℝ≥0} (ρ : ℝ)
    (hε : ε ≠ 0) (hρ : 0 ≤ ρ) :
    ∃ s : Finset M,
      Metric.IsSeparated ε (s : Set M) ∧
      (s : Set M) ⊆ Metric.closedBall o ρ ∧
      s.card = packingCount o ε ρ :=
  exists_optimal_blockCode o ε ρ hρ
    (hasFinitePacking_of_properSpace o hε)

/-- The packing theorem in a proper metric space, with the finiteness
    hypothesis discharged. The positive resolution field of the representation
    is exactly what licenses this. -/
theorem faithful_representation_addressable_of_properSpace
    [ProperSpace M] (o : M) (ε : ℝ≥0) (rep : RetainedRepresentation o ε)
    (β c hpack : ℝ)
    (hdemand : Tendsto (representedRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 hpack)) :
    Addressability.Addressable β c hpack :=
  faithful_representation_addressable o ε rep
    (hasFinitePacking_of_properSpace o rep.resolution_pos.ne') β c hpack
    hdemand hradial hpacking

/-- A zero-capacity host cannot carry positive retained growth under the
    hypotheses of the packing theorem. -/
theorem no_positive_growth_at_zero_capacity
    (o : M) (ε : ℝ≥0) (rep : RetainedRepresentation o ε)
    (hfinite : HasFinitePacking o ε)
    (β c : ℝ) (hβ : 0 < β)
    (hdemand : Tendsto (representedRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 0)) :
    False := by
  have hbound :=
    faithful_representation_addressable
      o ε rep hfinite β c 0 hdemand hradial hpacking
  unfold Addressability.Addressable at hbound
  nlinarith

end ActiveGeometry.Packing
