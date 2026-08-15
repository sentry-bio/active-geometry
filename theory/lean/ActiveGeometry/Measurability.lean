/-
  Growth-class gate identities (instrument, not a theory layer)
  ==============================================

  Elementary real identities behind the finite-sample growth-class test.
  See theory/MEASURABILITY.md.

  On a radial window of ratio r > 1, the unique exponential and the unique
  polynomial that match at both endpoints have log-occupancy difference

      f(t) = d log t - d (log r)/(r-1) (t-1),    t ∈ [1, r].

  This file proves:

  * the geometric-midpoint identity (Lemma 1.2);
  * the exact maximum-gap formula at the logarithmic mean (Proposition 1.3);
  * positivity of the span information Δ(r,d) for r > 1, d > 0.

  It does not formalize Hellinger distance, Le Cam's lemma, or the
  Poisson-increment testing model. Those remain paper proofs.
-/

import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

namespace ActiveGeometry.Measurability

open Real

/-- Log-occupancy difference of the endpoint-matched pair, on the scaled
    window t = ρ/ρ_min. -/
noncomputable def logGap (d r t : ℝ) : ℝ :=
  d * log t - d * (log r / (r - 1)) * (t - 1)

/-- Formal t-derivative of `logGap`. -/
noncomputable def logGapDeriv (d r t : ℝ) : ℝ :=
  d / t - d * (log r / (r - 1))

/-- Geometric-midpoint gap exponent (Lemma 1.2). -/
noncomputable def midpointExponent (d r : ℝ) : ℝ :=
  d * log r * (1 / 2 - 1 / (√r + 1))

/-- Exact maximum log-gap (Proposition 1.3). -/
noncomputable def spanInformation (d r : ℝ) : ℝ :=
  d * (log ((r - 1) / log r) - 1 + log r / (r - 1))

/-- Logarithmic mean of 1 and r: the unique critical point of `logGap`. -/
noncomputable def logMean (r : ℝ) : ℝ :=
  (r - 1) / log r

lemma logGap_one (d r : ℝ) (_hr : 1 < r) : logGap d r 1 = 0 := by
  unfold logGap
  simp [log_one]

lemma logGap_right (d r : ℝ) (hr : 1 < r) : logGap d r r = 0 := by
  have hr1 : r - 1 ≠ 0 := by linarith
  unfold logGap
  field_simp [hr1]
  ring

private lemma sqrt_factor (r : ℝ) (hr0 : 0 < r) :
    r - 1 = (√r - 1) * (√r + 1) := by
  nlinarith [sq_sqrt (le_of_lt hr0)]

private lemma midpoint_rate
    (r : ℝ) (hr : 1 < r) :
    log r / (r - 1) * (√r - 1) = log r / (√r + 1) := by
  have hr0 : 0 < r := lt_trans (by norm_num : (0 : ℝ) < 1) hr
  have hne : √r + 1 ≠ 0 := by
    have : 0 < √r := sqrt_pos.mpr hr0
    linarith
  have hsm : √r - 1 ≠ 0 := by
    intro h
    have hsqrt : √r = 1 := by linarith
    have : r = 1 := by
      calc
        r = (√r) ^ 2 := (sq_sqrt (le_of_lt hr0)).symm
        _ = 1 := by simp [hsqrt]
    linarith
  rw [sqrt_factor r hr0]
  field_simp [hne, hsm]

/-- Lemma 1.2: the midpoint identity. -/
theorem midpoint_exponent_eq (d r : ℝ) (hr : 1 < r) :
    logGap d r (√r) = midpointExponent d r := by
  have hr0 : 0 < r := lt_trans (by norm_num : (0 : ℝ) < 1) hr
  have hlog : log (√r) = log r / 2 := log_sqrt (le_of_lt hr0)
  have hquot := midpoint_rate r hr
  unfold logGap midpointExponent
  rw [hlog]
  calc
    d * (log r / 2) - d * (log r / (r - 1)) * (√r - 1)
      = d * (log r / 2) - d * (log r / (r - 1) * (√r - 1)) := by ring
    _ = d * (log r / 2) - d * (log r / (√r + 1)) := by rw [hquot]
    _ = d * log r * (1 / 2 - 1 / (√r + 1)) := by ring

/-- The formal derivative vanishes at the logarithmic mean. -/
theorem logGapDeriv_logMean (d r : ℝ) (hr : 1 < r) :
    logGapDeriv d r (logMean r) = 0 := by
  have hlog : 0 < log r := log_pos hr
  have hr1 : r - 1 ≠ 0 := by linarith
  unfold logGapDeriv logMean
  field_simp [hr1, ne_of_gt hlog]
  ring

private lemma logMean_shift
    (r : ℝ) (hr : 1 < r) :
    log r / (r - 1) * ((r - 1) / log r - 1) = 1 - log r / (r - 1) := by
  have hlog : 0 < log r := log_pos hr
  have hr1 : r - 1 ≠ 0 := by linarith
  field_simp [hr1, ne_of_gt hlog]

/-- Proposition 1.3: the gap at the logarithmic mean is the span information. -/
theorem spanInformation_eq_logGap_logMean (d r : ℝ) (hr : 1 < r) :
    logGap d r (logMean r) = spanInformation d r := by
  have ht := logMean_shift r hr
  unfold logGap spanInformation logMean
  calc
    d * log ((r - 1) / log r)
        - d * (log r / (r - 1)) * ((r - 1) / log r - 1)
      = d * log ((r - 1) / log r)
        - d * (log r / (r - 1) * ((r - 1) / log r - 1)) := by ring
    _ = d * log ((r - 1) / log r) - d * (1 - log r / (r - 1)) := by rw [ht]
    _ = d * (log ((r - 1) / log r) - 1 + log r / (r - 1)) := by ring

/-- Span information is positive for r > 1 and d > 0. -/
theorem spanInformation_pos (d r : ℝ) (hd : 0 < d) (hr : 1 < r) :
    0 < spanInformation d r := by
  have hlog : 0 < log r := log_pos hr
  have ht : 1 < logMean r := by
    unfold logMean
    exact (one_lt_div hlog).mpr (log_lt_sub_one_of_pos (by linarith) (ne_of_gt hr))
  have htpos : 0 < logMean r := lt_trans (by norm_num : (0 : ℝ) < 1) ht
  have hx : 0 < logMean r - 1 := sub_pos.mpr ht
  have hlogt : 2 * (logMean r - 1) / (logMean r + 1) < log (logMean r) := by
    have h := lt_log_one_add_of_pos hx
    have hden : logMean r - 1 + 2 = logMean r + 1 := by ring
    have harg : 1 + (logMean r - 1) = logMean r := by ring
    rw [hden, harg] at h
    exact h
  have hcmp :
      (logMean r - 1) / logMean r
        < 2 * (logMean r - 1) / (logMean r + 1) := by
    have hsum : 0 < logMean r + 1 := by linarith
    rw [div_lt_div_iff₀ htpos hsum]
    nlinarith [ht]
  have hmain : (logMean r - 1) / logMean r < log (logMean r) :=
    lt_trans hcmp hlogt
  have hrecip : log r / (r - 1) = (logMean r)⁻¹ := by
    unfold logMean
    field_simp [ne_of_gt hlog]
  have hform :
      logGap d r (logMean r)
        = d * (log (logMean r) - (logMean r - 1) / logMean r) := by
    unfold logGap
    rw [hrecip]
    field_simp [ne_of_gt htpos]
  have hdiff : 0 < log (logMean r) - (logMean r - 1) / logMean r :=
    sub_pos.mpr hmain
  have : 0 < d * (log (logMean r) - (logMean r - 1) / logMean r) :=
    mul_pos hd hdiff
  rw [← spanInformation_eq_logGap_logMean d r hr, hform]
  exact this

end ActiveGeometry.Measurability
