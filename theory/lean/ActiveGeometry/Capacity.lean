/-
  Capacity algebra of the addressability bound
  ============================================

  Pure real algebra downstream of

      Addressable β c hcap  :  β ≤ c · h_cap.

  Nothing here mentions a metric space, a packing number, DNA, or an alphabet.
  `Packing.lean` is what *proves* `Addressable` in the convergent-rate case.
  This file records what the inequality implies, and what it does not.

  The space-form formula `(n-1)√κ` is an algebraic identification of rates,
  written `hcap_eq_spaceForm`. It is not a theorem that the host is
  `ℍⁿ_κ`, and it is not the space-form classification (Theorem 7.1 of the
  spine is a paper sketch, not Lean).

  The equality case `β = c · h_cap` is a separate predicate
  `CapacitySaturated`. Combined with `hcap_eq_spaceForm` it becomes the
  face of the feasible region; that composition lives in `StateEquation.lean`.
-/

import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

namespace ActiveGeometry.Capacity

open Real

private lemma log2_pos : log 2 > 0 :=
  log_pos (by norm_num : (1 : ℝ) < 2)

/-- Convert a rate in bits per step to nats per step. -/
noncomputable def bitsToNats (h : ℝ) : ℝ := h * log 2

/-- The coordinate-free addressability condition `β ≤ c · h_cap`.
    Here `h_cap` is any independently established exponential host-capacity
    rate (fixed-resolution packing entropy, or volume entropy under extra
    hypotheses). -/
def Addressable (β c hcap : ℝ) : Prop := β ≤ c * hcap

/-- The equality case of addressability. This is not part of `Addressable`. -/
def CapacitySaturated (β c hcap : ℝ) : Prop := β = c * hcap

/-- Addressability efficiency `η = β / (c · h_cap)`. On the physical domain
    `β ≥ 0`, `c > 0`, `h_cap > 0`, the bound forces `η ≤ 1`. -/
noncomputable def efficiency (β c hcap : ℝ) : ℝ := β / (c * hcap)

/-- The real-hyperbolic space-form entropy `(n-1)√κ` (for `κ ≥ 0`).
    This is a formula for a rate, not a claim that a given host is `ℍⁿ_κ`. -/
noncomputable def spaceFormEntropy (n κval : ℝ) : ℝ :=
  (n - 1) * sqrt κval

/-- Algebraic identification of a host-capacity rate with the space-form
    value `(n-1)√κ`. This is not a classification theorem. -/
def hcap_eq_spaceForm (hcap n κval : ℝ) : Prop :=
  hcap = spaceFormEntropy n κval

/-- The least space-form curvature magnitude permitted by addressability on
    the positive physical domain. -/
noncomputable def curvatureFloor (β c n : ℝ) : ℝ :=
  (β / (c * (n - 1))) ^ 2

/-- Sectional-curvature magnitude in process-step units. This combination is
    invariant under a common rescaling of radial distance and curvature. -/
def normalizedCurvature (c κval : ℝ) : ℝ := c ^ 2 * κval

/-- The normalized curvature value on the saturation face. It is not raw
    sectional curvature unless the process-time gauge `c = 1` is chosen. -/
noncomputable def normalizedCurvatureAtSaturation (h n : ℝ) : ℝ :=
  (bitsToNats h / (n - 1)) ^ 2

/-- A positive bit rate gives a positive information-growth rate in nats. -/
theorem bitsToNats_pos (h : ℝ) (hh : 0 < h) :
    0 < bitsToNats h :=
  mul_pos hh log2_pos

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

/-- Saturation implies addressability; the converse is not assumed. -/
theorem saturated_is_addressable
    (β c hcap : ℝ) (hsaturated : CapacitySaturated β c hcap) :
    Addressable β c hcap := by
  unfold CapacitySaturated at hsaturated
  unfold Addressable
  exact hsaturated.le

/-- In a space-form chart, addressability imposes a curvature floor.
    This is the inequality form of the geometric state relation. -/
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
    addressability bound and a space-form identification of `h_cap`. -/
theorem addressable_spaceForm_floor
    (β c hcap n κval : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hbound : Addressable β c hcap)
    (hform : hcap_eq_spaceForm hcap n κval) :
    curvatureFloor β c n ≤ κval := by
  apply curvature_at_least_floor β c n κval hβ hc hn hκ
  unfold Addressable at hbound
  unfold hcap_eq_spaceForm at hform
  rw [hform] at hbound
  simpa [spaceFormEntropy, mul_assoc] using hbound

/-- Capacity saturation determines the raw sectional-curvature magnitude
    uniquely once the radial conversion `c` and ambient dimension `n` are
    fixed. -/
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

/-- Saturation fixes curvature only after the space-form identification has
    been supplied independently. -/
theorem saturated_spaceForm_eq_floor
    (β c hcap n κval : ℝ)
    (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hsaturated : CapacitySaturated β c hcap)
    (hform : hcap_eq_spaceForm hcap n κval) :
    κval = curvatureFloor β c n := by
  apply saturated_curvature_eq_floor β c n κval hc hn hκ
  unfold CapacitySaturated at hsaturated
  unfold hcap_eq_spaceForm at hform
  rw [hform] at hsaturated
  simpa [spaceFormEntropy, mul_assoc] using hsaturated

/-- The curvature floor itself saturates space-form capacity. -/
theorem floor_saturates_capacity
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    c * (n - 1) * sqrt (curvatureFloor β c n) = β := by
  unfold curvatureFloor
  have hden : 0 < c * (n - 1) := mul_pos hc (by linarith)
  have hn0 : n - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_gt hn)
  rw [sqrt_sq (div_pos hβ hden).le]
  field_simp [hden.ne', hn0]

/-- Multiplying the raw curvature-magnitude floor by `c²` removes radial-unit
    dependence and recovers the normalized equality-case expression. -/
theorem normalized_floor_eq_saturation
    (h c n : ℝ) (hc : c ≠ 0) (hn : n ≠ 1) :
    normalizedCurvature c (curvatureFloor (bitsToNats h) c n) =
      normalizedCurvatureAtSaturation h n := by
  unfold normalizedCurvature curvatureFloor normalizedCurvatureAtSaturation
  have hn0 : n - 1 ≠ 0 := sub_ne_zero.mpr hn
  field_simp [hc, hn0]

/-- The process-time gauge `c = 1` turns the raw curvature-magnitude floor
    into the familiar formula. -/
theorem process_time_gauge
    (h n : ℝ) :
    curvatureFloor (bitsToNats h) 1 n =
      normalizedCurvatureAtSaturation h n := by
  simp [curvatureFloor, normalizedCurvatureAtSaturation]

/-- Normalized curvature magnitude is invariant when radial distance is
    rescaled by a nonzero factor `a` and raw curvature magnitude by `a⁻²`. -/
theorem normalized_curvature_scale_invariant
    (a c κval : ℝ) (ha : a ≠ 0) :
    normalizedCurvature (a * c) (κval / a ^ 2) =
      normalizedCurvature c κval := by
  unfold normalizedCurvature
  field_simp [ha]

/-- At ambient dimension 2, the equality-case formula is `(h ln 2)²`. -/
theorem normalizedCurvatureAtSaturation_two (h : ℝ) :
    normalizedCurvatureAtSaturation h 2 = (bitsToNats h) ^ 2 := by
  simp [normalizedCurvatureAtSaturation]; ring

theorem normalizedCurvatureAtSaturation_two_pos (h : ℝ) (hh : 0 < h) :
    0 < normalizedCurvatureAtSaturation h 2 := by
  rw [normalizedCurvatureAtSaturation_two]
  exact sq_pos_of_pos (bitsToNats_pos h hh)

/-- The equality-case formula is strictly increasing in the bit rate. -/
theorem normalizedCurvatureAtSaturation_mono_h
    (h₁ h₂ : ℝ) (h1pos : 0 < h₁) (hlt : h₁ < h₂) :
    normalizedCurvatureAtSaturation h₁ 2 <
      normalizedCurvatureAtSaturation h₂ 2 := by
  simp only [normalizedCurvatureAtSaturation_two, bitsToNats]
  exact sq_lt_sq' (by nlinarith [log2_pos]) (by nlinarith [log2_pos])

/-- The equality-case formula is strictly decreasing in ambient dimension.
    This is algebraic monotonicity, not a selection principle for `n = 2`. -/
theorem normalizedCurvatureAtSaturation_anti_n
    (h n₁ n₂ : ℝ) (hh : 0 < h)
    (hn1 : 1 < n₁) (hn2 : 1 < n₂) (hlt : n₁ < n₂) :
    normalizedCurvatureAtSaturation h n₂ <
      normalizedCurvatureAtSaturation h n₁ := by
  unfold normalizedCurvatureAtSaturation bitsToNats
  have hnum := mul_pos hh log2_pos
  have hd1 : n₁ - 1 > 0 := sub_pos.mpr hn1
  have hd2 : n₂ - 1 > 0 := sub_pos.mpr hn2
  exact sq_lt_sq' (by nlinarith [div_pos hnum hd2, div_pos hnum hd1])
    (div_lt_div_of_pos_left hnum hd1 (by linarith))

theorem normalizedCurvatureAtSaturation_lt_of_dim_gt_two
    (h n : ℝ) (hh : 0 < h) (hn : 2 < n) :
    normalizedCurvatureAtSaturation h n <
      normalizedCurvatureAtSaturation h 2 :=
  normalizedCurvatureAtSaturation_anti_n
    h 2 n hh (by norm_num) (by linarith) hn

/-- Scaling the bit rate scales the equality-case curvature by the square. -/
theorem normalizedCurvatureAtSaturation_mul_h (h a : ℝ) :
    normalizedCurvatureAtSaturation (a * h) 2 =
      a ^ 2 * normalizedCurvatureAtSaturation h 2 := by
  simp [normalizedCurvatureAtSaturation_two, bitsToNats]; ring

/-- At `n = 2`, radius times space-form entropy of the equality-case
    curvature recovers radius times the information rate. -/
theorem radial_infoRate_at_saturation_two (h r : ℝ) (hh : 0 < h) :
    r * sqrt (normalizedCurvatureAtSaturation h 2) = r * bitsToNats h := by
  rw [normalizedCurvatureAtSaturation_two,
    sqrt_sq (bitsToNats_pos h hh).le]

end ActiveGeometry.Capacity
