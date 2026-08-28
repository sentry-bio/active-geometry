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

  Cite `Packing.addressability_limit` for the bound. Cite
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
    normalizedCurvature c κval = idealNormalizedCurvature h n := by
  have hkappa : κval = curvatureFloor (bitsToNats h) c n :=
    saturated_spaceForm_eq_floor
      (bitsToNats h) c hcap n κval hc hn hκ hsaturated hform
  rw [hkappa]
  exact normalized_floor_eq_ideal h c n hc.ne' (ne_of_gt hn)

/-- Signed gap between an information rate and space-form entropy.
    Vanishes on the equality-case face (at process-time gauge `c = 1`). -/
noncomputable def rateMismatch (h κval n : ℝ) : ℝ :=
  bitsToNats h - spaceFormEntropy n κval

/-- Squared rate gap. Non-negative by construction. -/
noncomputable def rateMismatchSq (h κval n : ℝ) : ℝ :=
  (rateMismatch h κval n) ^ 2

/-- Squared gap of square roots of curvature magnitudes. -/
noncomputable def sqrtCurvatureMismatch (κval κstar : ℝ) : ℝ :=
  (sqrt κval - sqrt κstar) ^ 2

theorem rateMismatchSq_nonneg (h κval n : ℝ) : 0 ≤ rateMismatchSq h κval n :=
  sq_nonneg _

theorem sqrtCurvatureMismatch_nonneg (κval κstar : ℝ) :
    0 ≤ sqrtCurvatureMismatch κval κstar :=
  sq_nonneg _

private lemma sqrt_idealNormalizedCurvature
    (h n : ℝ) (hh : 0 < h) (hn : 1 < n) :
    sqrt (idealNormalizedCurvature h n) = bitsToNats h / (n - 1) := by
  unfold idealNormalizedCurvature
  exact sqrt_sq (div_pos (bitsToNats_pos h hh) (by linarith)).le

/-- On positive curvature magnitudes, the square-root mismatch vanishes
    exactly at equality. -/
theorem sqrtCurvatureMismatch_zero_iff (κval κstar : ℝ)
    (hκpos : 0 < κval) (hκspos : 0 < κstar) :
    sqrtCurvatureMismatch κval κstar = 0 ↔ κval = κstar := by
  unfold sqrtCurvatureMismatch
  constructor
  · intro hV
    have := sq_eq_zero_iff.mp hV
    calc κval = (sqrt κval) ^ 2 := (sq_sqrt hκpos.le).symm
      _ = (sqrt κstar) ^ 2 := by rw [show sqrt κval = sqrt κstar by linarith]
      _ = κstar := sq_sqrt hκspos.le
  · intro heq; rw [heq, sub_self]; ring

/-- The squared rate gap and the square-root curvature mismatch are the same
    object up to the dimensional factor `(n-1)²`. -/
theorem rateMismatchSq_eq_scaled_sqrtMismatch (h κval n : ℝ)
    (hh : 0 < h) (_hκ : 0 < κval) (hn : 1 < n) :
    rateMismatchSq h κval n =
      (n - 1) ^ 2 * sqrtCurvatureMismatch κval (idealNormalizedCurvature h n) := by
  have hn1 : n - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_gt hn)
  unfold rateMismatchSq rateMismatch sqrtCurvatureMismatch spaceFormEntropy
  rw [sqrt_idealNormalizedCurvature h n hh hn]
  have hscale :
      (n - 1) * (sqrt κval - bitsToNats h / (n - 1)) =
        (n - 1) * sqrt κval - bitsToNats h := by
    field_simp [hn1]
  calc
    (bitsToNats h - (n - 1) * sqrt κval) ^ 2
        = ((n - 1) * sqrt κval - bitsToNats h) ^ 2 := by ring
    _ = ((n - 1) * (sqrt κval - bitsToNats h / (n - 1))) ^ 2 := by rw [hscale]
    _ = (n - 1) ^ 2 * (sqrt κval - bitsToNats h / (n - 1)) ^ 2 := by ring

private lemma idealNormalizedCurvature_pos
    (h n : ℝ) (hh : 0 < h) (hn : 1 < n) :
    0 < idealNormalizedCurvature h n :=
  pow_pos (div_pos (bitsToNats_pos h hh) (by linarith)) 2

/-- The squared rate gap vanishes exactly at the equality-case curvature.
    A corollary of the scaling identity, not an independent computation. -/
theorem rateMismatchSq_zero_iff (h κval n : ℝ)
    (hh : 0 < h) (hκpos : 0 < κval) (hn : 1 < n) :
    rateMismatchSq h κval n = 0 ↔ κval = idealNormalizedCurvature h n := by
  have hcrit_pos := idealNormalizedCurvature_pos h n hh hn
  have hfac : (n - 1) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.mpr (ne_of_gt hn))
  rw [rateMismatchSq_eq_scaled_sqrtMismatch h κval n hh hκpos hn, mul_eq_zero,
    sqrtCurvatureMismatch_zero_iff κval (idealNormalizedCurvature h n)
      hκpos hcrit_pos]
  simp [hfac]

/-- At process-time gauge `c = 1`, the space-form entropy of the equality-case
    curvature recovers the information rate. This is not a statement about
    the derivative of a physical potential. -/
theorem rateMismatch_zero_at_face
    (h n : ℝ) (hh : 0 < h) (hn : 1 < n) :
    rateMismatch h (idealNormalizedCurvature h n) n = 0 := by
  have hsat :=
    floor_saturates_capacity (bitsToNats h) 1 n (bitsToNats_pos h hh)
      (by norm_num) hn
  rw [process_time_gauge] at hsat
  unfold rateMismatch spaceFormEntropy
  linarith

end ActiveGeometry.StateEquation
