/-
  The equality-case face (optional)
  =================================

  One inequality defines a feasible region: `Addressable β c h_cap`.
  This file is the *face* of that region in a space-form chart — the
  composition of `CapacitySaturated` with `hcap_eq_spaceForm` — plus
  algebraic diagnostics of the gap to that face.

  It is not a second primitive, not an evolution law, and not a Lyapunov
  theorem. The functions `rateMismatchSq` and `sqrtCurvatureMismatch` are
  positive-definite gap diagnostics; their unique-zero results do not prove
  stability without an explicit dynamics and a proof that the diagnostic
  decreases along trajectories.

  Cite `Packing.convergent_rate_addressability_limit` for the Lean-checked
  ordinary-limit bound. Cite
  `normalized_state_equation` only when saturation and a space-form
  identification have been supplied independently.
-/

import ActiveGeometry.Capacity
import Mathlib.Tactic

namespace ActiveGeometry.StateEquation

open Real
open Capacity

/-- The state equation is the conditional equality case of the kernel:
    capacity saturation plus a space-form identification of `h_cap`. -/
theorem normalized_state_equation
    (h c hcap n κval : ℝ)
    (hc : 0 < c) (hn : 1 < n) (hκ : 0 ≤ κval)
    (hsaturated : CapacitySaturated (bitsToNats h) c hcap)
    (hform : hcap_eq_spaceForm hcap n κval) :
    normalizedCurvature c κval = normalizedCurvatureAtSaturation h n := by
  have hkappa : κval = curvatureFloor (bitsToNats h) c n :=
    saturated_spaceForm_eq_floor
      (bitsToNats h) c hcap n κval hc hn hκ hsaturated hform
  rw [hkappa]
  exact normalized_floor_eq_saturation h c n hc.ne' (ne_of_gt hn)

/-- Signed gap between demand and space-form capacity.
    It vanishes exactly on the equality-case face
    `β = c · (n-1)√κ`. No gauge has been fixed. -/
noncomputable def rateMismatch (β c κval n : ℝ) : ℝ :=
  β - c * spaceFormEntropy n κval

/-- Squared rate gap. Non-negative by construction. -/
noncomputable def rateMismatchSq (β c κval n : ℝ) : ℝ :=
  (rateMismatch β c κval n) ^ 2

/-- Squared gap of square roots of curvature magnitudes. -/
noncomputable def sqrtCurvatureMismatch (κval κstar : ℝ) : ℝ :=
  (sqrt κval - sqrt κstar) ^ 2

theorem rateMismatchSq_nonneg (β c κval n : ℝ) :
    0 ≤ rateMismatchSq β c κval n :=
  sq_nonneg _

theorem sqrtCurvatureMismatch_nonneg (κval κstar : ℝ) :
    0 ≤ sqrtCurvatureMismatch κval κstar :=
  sq_nonneg _

private lemma sqrt_curvatureFloor
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    sqrt (curvatureFloor β c n) = β / (c * (n - 1)) := by
  unfold curvatureFloor
  exact sqrt_sq (div_pos hβ (mul_pos hc (by linarith))).le

/-- On non-negative curvature magnitudes, the square-root mismatch vanishes
    exactly at equality. -/
theorem sqrtCurvatureMismatch_zero_iff (κval κstar : ℝ)
    (hκ : 0 ≤ κval) (hκstar : 0 ≤ κstar) :
    sqrtCurvatureMismatch κval κstar = 0 ↔ κval = κstar := by
  unfold sqrtCurvatureMismatch
  constructor
  · intro hV
    have := sq_eq_zero_iff.mp hV
    calc κval = (sqrt κval) ^ 2 := (sq_sqrt hκ).symm
      _ = (sqrt κstar) ^ 2 := by rw [show sqrt κval = sqrt κstar by linarith]
      _ = κstar := sq_sqrt hκstar
  · intro heq; rw [heq, sub_self]; ring

/-- The squared rate gap and the square-root curvature mismatch are the same
    object up to the dimensional factor `(n-1)²`. -/
theorem rateMismatchSq_eq_scaled_sqrtMismatch (β c κval n : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    rateMismatchSq β c κval n =
      (c * (n - 1)) ^ 2 *
        sqrtCurvatureMismatch κval (curvatureFloor β c n) := by
  have hscale : c * (n - 1) ≠ 0 :=
    (mul_pos hc (by linarith)).ne'
  unfold rateMismatchSq rateMismatch sqrtCurvatureMismatch spaceFormEntropy
  rw [sqrt_curvatureFloor β c n hβ hc hn]
  have hfactor :
      (c * (n - 1)) * (sqrt κval - β / (c * (n - 1))) =
        c * (n - 1) * sqrt κval - β := by
    rw [mul_sub, mul_div_cancel₀ β hscale]
  calc
    (β - c * ((n - 1) * sqrt κval)) ^ 2
        = (c * (n - 1) * sqrt κval - β) ^ 2 := by ring
    _ = ((c * (n - 1)) *
          (sqrt κval - β / (c * (n - 1)))) ^ 2 := by rw [hfactor]
    _ = (c * (n - 1)) ^ 2 *
          (sqrt κval - β / (c * (n - 1))) ^ 2 := by ring

private lemma curvatureFloor_pos
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    0 < curvatureFloor β c n :=
  pow_pos (div_pos hβ (mul_pos hc (by linarith))) 2

/-- The squared rate gap vanishes exactly at the equality-case curvature.
    A corollary of the scaling identity, not an independent computation. -/
theorem rateMismatchSq_zero_iff (β c κval n : ℝ)
    (hβ : 0 < β) (hc : 0 < c) (hκ : 0 ≤ κval) (hn : 1 < n) :
    rateMismatchSq β c κval n = 0 ↔ κval = curvatureFloor β c n := by
  have hfloor_pos := curvatureFloor_pos β c n hβ hc hn
  have hfac : (c * (n - 1)) ^ 2 ≠ 0 :=
    pow_ne_zero 2 (mul_ne_zero hc.ne' (sub_ne_zero.mpr (ne_of_gt hn)))
  rw [rateMismatchSq_eq_scaled_sqrtMismatch β c κval n hβ hc hn, mul_eq_zero,
    sqrtCurvatureMismatch_zero_iff κval (curvatureFloor β c n)
      hκ hfloor_pos.le]
  simp [hfac]

/-- Space-form capacity at the curvature floor equals demand.
    This is not a statement about the derivative of a physical potential. -/
theorem rateMismatch_zero_at_floor
    (β c n : ℝ) (hβ : 0 < β) (hc : 0 < c) (hn : 1 < n) :
    rateMismatch β c (curvatureFloor β c n) n = 0 := by
  have hsat :=
    floor_saturates_capacity β c n hβ hc hn
  unfold rateMismatch spaceFormEntropy
  linarith

end ActiveGeometry.StateEquation
