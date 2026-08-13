/-
  The Addressability Limit
  ========================

  Algebraic consequences of the coordinate-free packing bound

      β ≤ c · h_vol

  where β is retained-information growth (nats / generative step), c is
  radial distance / generative step, and h_vol is host volume entropy
  (nats / radial distance).

  For β > 0, c > 0, n > 1, and κ ≥ 0, an isotropic hyperbolic host has
  h_vol = (n - 1) √κ, so the bound gives

      κ ≥ (β / (c (n - 1)))².

  Equality is the capacity-saturating case. For β = h ln 2 and normalized
  curvature κ̄ = c²κ, this becomes

      κ̄ = (h ln 2 / (n - 1))².

  This file formalizes the algebra after the packing bound is supplied as a
  hypothesis. It does not formalize the metric packing theorem, the space-form
  classification, or a dynamics toward saturation.
-/

import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

namespace ActiveGeometry.Addressability

open Real

private lemma log2_pos : log 2 > 0 :=
  log_pos (by norm_num : (1 : ℝ) < 2)

/-- Convert a transmitted rate in bits per step to nats per step. -/
noncomputable def bitsToNats (h : ℝ) : ℝ := h * log 2

/-- The coordinate-free addressability condition β ≤ c h_vol. -/
def Addressable (β c hvol : ℝ) : Prop := β ≤ c * hvol

/-- Addressability efficiency. Its physical domain has β ≥ 0, c > 0,
    h_vol > 0; under the addressability bound it is at most one. -/
noncomputable def efficiency (β c hvol : ℝ) : ℝ := β / (c * hvol)

/-- Isotropic hyperbolic volume entropy (for κ ≥ 0). -/
noncomputable def isotropicEntropy (n κval : ℝ) : ℝ :=
  (n - 1) * sqrt κval

/-- The least isotropic sectional-curvature magnitude permitted by
    addressability on the positive physical domain. -/
noncomputable def curvatureFloor (β c n : ℝ) : ℝ :=
  (β / (c * (n - 1))) ^ 2

/-- Sectional-curvature magnitude in process-step units. This combination is
    invariant under a common rescaling of radial distance and curvature. -/
def normalizedCurvature (c κval : ℝ) : ℝ := c ^ 2 * κval

/-- The familiar state-equation expression is normalized curvature, not raw
    sectional curvature unless the process-time gauge c = 1 is chosen. -/
noncomputable def idealNormalizedCurvature (h n : ℝ) : ℝ :=
  (bitsToNats h / (n - 1)) ^ 2

/-- Positive retained-information growth at finite positive radial rate forces
    positive host entropy. -/
theorem addressability_forces_positive_entropy
    (β c hvol : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hbound : Addressable β c hvol) :
    0 < hvol := by
  unfold Addressable at hbound
  by_contra hnot
  have hh : hvol ≤ 0 := le_of_not_gt hnot
  have hprod : c * hvol ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (le_of_lt hc) hh
  linarith

/-- The addressability efficiency cannot exceed one. -/
theorem efficiency_le_one
    (β c hvol : ℝ)
    (hc : 0 < c) (hh : 0 < hvol) (hbound : Addressable β c hvol) :
    efficiency β c hvol ≤ 1 := by
  unfold efficiency
  rw [div_le_iff₀ (mul_pos hc hh)]
  simpa [Addressable] using hbound

/-- In an isotropic hyperbolic model, addressability imposes a curvature
    floor. This is the inequality form of the geometric state relation. -/
theorem curvature_at_least_floor
    (β c n κval : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hcapacity : β ≤ c * (n - 1) * sqrt κval) :
    curvatureFloor β c n ≤ κval := by
  unfold curvatureFloor
  have hden : 0 < c * (n - 1) := mul_pos hc (by linarith)
  have hquot : β / (c * (n - 1)) ≤ sqrt κval := by
    rw [div_le_iff₀ hden]
    calc
      β ≤ c * (n - 1) * sqrt κval := hcapacity
      _ = sqrt κval * (c * (n - 1)) := by ring
  have hquot_nonneg : 0 ≤ β / (c * (n - 1)) :=
    le_of_lt (div_pos hβ hden)
  have hsqrt_nonneg : 0 ≤ sqrt κval := sqrt_nonneg κval
  have hsqrt_sq : (sqrt κval) ^ 2 = κval := sq_sqrt hκ
  nlinarith

/-- Capacity saturation determines the raw sectional-curvature magnitude
    uniquely once the radial conversion c and ambient dimension n are fixed. -/
theorem saturated_curvature_eq_floor
    (β c n κval : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hsaturated : β = c * (n - 1) * sqrt κval) :
    κval = curvatureFloor β c n := by
  have hden : 0 < c * (n - 1) := mul_pos hc (by linarith)
  have hsqrt : sqrt κval = β / (c * (n - 1)) := by
    rw [eq_div_iff hden.ne']
    calc
      sqrt κval * (c * (n - 1))
          = c * (n - 1) * sqrt κval := by ring
      _ = β := hsaturated.symm
  calc
    κval = (sqrt κval) ^ 2 := (sq_sqrt hκ).symm
    _ = (β / (c * (n - 1))) ^ 2 := by rw [hsqrt]
    _ = curvatureFloor β c n := rfl

/-- The curvature floor itself saturates isotropic capacity. -/
theorem floor_saturates_capacity
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    c * (n - 1) * sqrt (curvatureFloor β c n) = β := by
  unfold curvatureFloor
  have hden : 0 < c * (n - 1) := mul_pos hc (by linarith)
  rw [sqrt_sq (div_pos hβ hden).le]
  field_simp [hden.ne']

/-- Multiplying the raw curvature-magnitude floor by c² removes radial-unit
    dependence and recovers the normalized state-equation expression. -/
theorem normalized_floor_eq_ideal
    (h c n : ℝ) (hc : c ≠ 0) (hn : n ≠ 1) :
    normalizedCurvature c (curvatureFloor (bitsToNats h) c n) =
      idealNormalizedCurvature h n := by
  unfold normalizedCurvature curvatureFloor idealNormalizedCurvature
  have hn0 : n - 1 ≠ 0 := sub_ne_zero.mpr hn
  field_simp [hc, hn0]
  <;> ring

/-- The process-time gauge c = 1 turns the raw curvature-magnitude floor into
    the familiar formula. -/
theorem process_time_gauge
    (h n : ℝ) :
    curvatureFloor (bitsToNats h) 1 n =
      idealNormalizedCurvature h n := by
  simp [curvatureFloor, idealNormalizedCurvature]

/-- Normalized curvature magnitude is invariant when radial distance is
    rescaled by a nonzero factor a and raw curvature magnitude by a⁻². -/
theorem normalized_curvature_scale_invariant
    (a c κval : ℝ) (ha : a ≠ 0) :
    normalizedCurvature (a * c) (κval / a ^ 2) =
      normalizedCurvature c κval := by
  unfold normalizedCurvature
  field_simp [ha]
  <;> ring

/-- A positive bit rate gives a positive information-growth rate in nats. -/
theorem bitsToNats_pos (h : ℝ) (hh : 0 < h) :
    0 < bitsToNats h :=
  mul_pos hh log2_pos

end ActiveGeometry.Addressability
