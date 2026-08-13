/-
  The Addressability Limit
  ========================

  Algebraic consequences of the coordinate-free capacity bound

      β ≤ c · h_cap

  where β is retained-information growth (nats / generative step), c is
  radial distance / generative step, and h_cap is independently established
  host packing entropy (or volume entropy under additional hypotheses), in
  nats / radial distance.

  For β > 0, c > 0, n > 1, and κ ≥ 0, an isotropic hyperbolic host has
  h_cap = h_vol = (n - 1) √κ, so the bound gives

      κ ≥ (β / (c (n - 1)))².

  Capacity saturation and isotropic realization are separate hypotheses.
  When both hold, β = h ln 2, and normalized curvature κ̄ = c²κ, the
  equality becomes

      κ̄ = (h ln 2 / (n - 1))².

  The predicates Addressable, CapacitySaturated, and IsotropicHyperbolic keep
  these levels explicit. This file formalizes the algebra after the capacity
  bound is supplied as a hypothesis. `Packing.lean` proves that hypothesis from
  metric packing in the convergent-rate case. Neither file formalizes the
  space-form classification or a dynamics toward saturation.
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

/-- The coordinate-free addressability condition β ≤ c hcap. Here hcap is
    any independently established exponential host-capacity rate, such as
    fixed-resolution packing entropy or, under additional hypotheses, volume
    entropy. -/
def Addressable (β c hcap : ℝ) : Prop := β ≤ c * hcap

/-- Capacity saturation is a separate condition, not part of addressability. -/
def CapacitySaturated (β c hcap : ℝ) : Prop := β = c * hcap

/-- Addressability efficiency. Its physical domain has β ≥ 0, c > 0,
    h_cap > 0; under the addressability bound it is at most one. -/
noncomputable def efficiency (β c hcap : ℝ) : ℝ := β / (c * hcap)

/-- Isotropic hyperbolic volume entropy (for κ ≥ 0). -/
noncomputable def isotropicEntropy (n κval : ℝ) : ℝ :=
  (n - 1) * sqrt κval

/-- An isotropic hyperbolic realization identifies independently specified
    host entropy with the space-form value (n - 1) √κ. -/
def IsotropicHyperbolic (hcap n κval : ℝ) : Prop :=
  hcap = isotropicEntropy n κval

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
    (β c hcap : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hbound : Addressable β c hcap) :
    0 < hcap := by
  unfold Addressable at hbound
  by_contra hnot
  have hh : hcap ≤ 0 := le_of_not_gt hnot
  have hprod : c * hcap ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (le_of_lt hc) hh
  linarith

/-- The addressability efficiency cannot exceed one. -/
theorem efficiency_le_one
    (β c hcap : ℝ)
    (hc : 0 < c) (hh : 0 < hcap) (hbound : Addressable β c hcap) :
    efficiency β c hcap ≤ 1 := by
  unfold efficiency
  rw [div_le_iff₀ (mul_pos hc hh)]
  simpa [Addressable] using hbound

/-- Saturation implies addressability, but the converse is not assumed. -/
theorem saturated_is_addressable
    (β c hcap : ℝ) (hsaturated : CapacitySaturated β c hcap) :
    Addressable β c hcap := by
  unfold CapacitySaturated at hsaturated
  unfold Addressable
  exact hsaturated.le

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

/-- Composition of the two independent ingredients: the coordinate-free
    addressability bound and an isotropic hyperbolic realization. -/
theorem isotropic_curvature_at_least_floor
    (β c hcap n κval : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hbound : Addressable β c hcap)
    (hisotropic : IsotropicHyperbolic hcap n κval) :
    curvatureFloor β c n ≤ κval := by
  apply curvature_at_least_floor β c n κval hβ hc hn hκ
  unfold Addressable at hbound
  unfold IsotropicHyperbolic at hisotropic
  rw [hisotropic] at hbound
  simpa [isotropicEntropy, mul_assoc] using hbound

/-- Capacity saturation determines the raw sectional-curvature magnitude
    uniquely once the radial conversion c and ambient dimension n are fixed. -/
theorem saturated_curvature_eq_floor
    (β c n κval : ℝ)
    (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
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

/-- Saturation fixes curvature only after the isotropic realization has been
    supplied independently. -/
theorem isotropic_saturation_curvature_eq_floor
    (β c hcap n κval : ℝ)
    (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hsaturated : CapacitySaturated β c hcap)
    (hisotropic : IsotropicHyperbolic hcap n κval) :
    κval = curvatureFloor β c n := by
  apply saturated_curvature_eq_floor β c n κval hc hn hκ
  unfold CapacitySaturated at hsaturated
  unfold IsotropicHyperbolic at hisotropic
  rw [hisotropic] at hsaturated
  simpa [isotropicEntropy, mul_assoc] using hsaturated

/-- The curvature floor itself saturates isotropic capacity. -/
theorem floor_saturates_capacity
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    c * (n - 1) * sqrt (curvatureFloor β c n) = β := by
  unfold curvatureFloor
  have hden : 0 < c * (n - 1) := mul_pos hc (by linarith)
  have hn0 : n - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_gt hn)
  rw [sqrt_sq (div_pos hβ hden).le]
  field_simp [hden.ne', hn0]

/-- Multiplying the raw curvature-magnitude floor by c² removes radial-unit
    dependence and recovers the normalized state-equation expression. -/
theorem normalized_floor_eq_ideal
    (h c n : ℝ) (hc : c ≠ 0) (hn : n ≠ 1) :
    normalizedCurvature c (curvatureFloor (bitsToNats h) c n) =
      idealNormalizedCurvature h n := by
  unfold normalizedCurvature curvatureFloor idealNormalizedCurvature
  have hn0 : n - 1 ≠ 0 := sub_ne_zero.mpr hn
  field_simp [hc, hn0]

/-- The state equation is the conditional equality case of the kernel:
    capacity saturation plus an isotropic hyperbolic realization. -/
theorem normalized_state_equation
    (h c hcap n κval : ℝ)
    (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hsaturated : CapacitySaturated (bitsToNats h) c hcap)
    (hisotropic : IsotropicHyperbolic hcap n κval) :
    normalizedCurvature c κval = idealNormalizedCurvature h n := by
  have hkappa : κval = curvatureFloor (bitsToNats h) c n :=
    isotropic_saturation_curvature_eq_floor
      (bitsToNats h) c hcap n κval hc hn hκ hsaturated hisotropic
  rw [hkappa]
  exact normalized_floor_eq_ideal h c n hc.ne' (ne_of_gt hn)

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

/-- A positive bit rate gives a positive information-growth rate in nats. -/
theorem bitsToNats_pos (h : ℝ) (hh : 0 < h) :
    0 < bitsToNats h :=
  mul_pos hh log2_pos

end ActiveGeometry.Addressability
