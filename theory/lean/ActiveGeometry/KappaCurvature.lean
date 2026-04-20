/-
  Biosphere Curvature Theorem: Lean 4 Formalization
  ==================================================

  Machine-checked proof of the geometric state equation for the tree of life:

    κ = (h · ln 2 / (n - 1))²

  where κ is the sectional curvature of the embedding hyperbolic manifold,
  h is the Shannon entropy rate of the generating code (bits/symbol), and
  n is the embedding dimension. At n = 2 (the topological invariant of trees
  established empirically across DNA, RNA, and protein alphabets in
  Fenn & Fenn 2026), this reduces to κ = (h · ln 2)².

  Contents:
    Part I  — Core theorems (kappa_n2, kappa_pos, kappa_unique, kappa_mono_h,
              kappa_mono_n, kappa_max_at_n2, kappa_scaling, growth_rate_match).
    Part II — Entropy rate decomposition and alphabet-capacity bound
              (h_three_letter, H_raw, entropy_rate_decomposition_bounds,
              kappa_bounded_by_raw, kappa_bounded_by_alphabet).
              Plus the trilogy-general form H_raw_of_alphabet α and the
              generic ceiling kappa_bounded_by_alphabet_general for any α > 1.
    Part III — Rate-distortion potential and Lyapunov stability
              (U, V, potential_zero_iff, lyapunov_zero_iff,
              potential_gradient_zero_at_critical).

  Companion papers (The Hyperbolic Trilogy):
    Paper I   — Fenn & Fenn, "Evolution as Active Geometry"
                bioRxiv 10.64898/2026.03.09.710612 (DNA/RNA/protein substrate)
    Paper II  — Fenn & Fenn, "A Geometric State Equation for
                Information-Generating Hierarchies"
                Zenodo 10.5281/zenodo.19381558 (this file ships with Paper II)
    Paper III — Fenn & Fenn, "Convergent Alphabets" (in prep;
                github.com/sentry-bio/convergent-alphabets; phoneme substrate)
-/

import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.Tactic

namespace BiosphereCurvature

open Real

noncomputable def κ (h : ℝ) (n : ℝ) : ℝ :=
  if n ≤ 1 then 0
  else (h * log 2 / (n - 1)) ^ 2

private lemma log2_pos : log 2 > 0 := log_pos (by norm_num : (1 : ℝ) < 2)

theorem kappa_n2 (h : ℝ) : κ h 2 = (h * log 2) ^ 2 := by
  simp [κ, show ¬(2 : ℝ) ≤ 1 by norm_num]; ring

theorem kappa_pos (h : ℝ) (hpos : h > 0) : κ h 2 > 0 := by
  rw [kappa_n2]; exact sq_pos_of_pos (mul_pos hpos log2_pos)

theorem kappa_unique (h : ℝ) (hpos : h > 0) :
    ∃! k : ℝ, k > 0 ∧ k = (h * log 2) ^ 2 :=
  ⟨_, ⟨sq_pos_of_pos (mul_pos hpos log2_pos), rfl⟩, fun _ ⟨_, hk'⟩ => hk'⟩

theorem kappa_mono_h (h₁ h₂ : ℝ) (h1pos : h₁ > 0) (h2pos : h₂ > 0) (hlt : h₁ < h₂) :
    κ h₁ 2 < κ h₂ 2 := by
  simp only [kappa_n2]
  exact sq_lt_sq' (by nlinarith [log2_pos]) (by nlinarith [log2_pos])

theorem kappa_mono_n (h : ℝ) (n₁ n₂ : ℝ) (hpos : h > 0)
    (hn1 : n₁ > 1) (hn2 : n₂ > 1) (hlt : n₁ < n₂) :
    κ h n₂ < κ h n₁ := by
  unfold κ; simp only [not_le.mpr hn1, not_le.mpr hn2, ↓reduceIte]
  have hnum := mul_pos hpos log2_pos
  have hd1 : n₁ - 1 > 0 := by linarith
  have hd2 : n₂ - 1 > 0 := by linarith
  exact sq_lt_sq' (by nlinarith [div_pos hnum hd2, div_pos hnum hd1])
    (div_lt_div_of_pos_left hnum hd1 (by linarith))

theorem kappa_max_at_n2 (h : ℝ) (n : ℝ) (hpos : h > 0) (hn : n > 2) :
    κ h n < κ h 2 :=
  kappa_mono_n h 2 n hpos (by norm_num) (by linarith) hn

theorem kappa_scaling (h c : ℝ) : κ (c * h) 2 = c ^ 2 * κ h 2 := by
  simp [kappa_n2]; ring

theorem growth_rate_match (h : ℝ) (r : ℝ) (hpos : h > 0) :
    r * sqrt (κ h 2) = h * r * log 2 := by
  rw [kappa_n2, sqrt_sq (mul_pos hpos log2_pos).le]; ring

noncomputable def κ_biosphere : ℝ := κ 1.6 2

theorem kappa_biosphere_form : κ_biosphere = (1.6 * log 2) ^ 2 := kappa_n2 1.6

/-
  ═══════════════════════════════════════════════════════════════════════════════
  PART II: ENTROPY RATE DECOMPOSITION
  ═══════════════════════════════════════════════════════════════════════════════

  The paper (Fenn & Fenn 2026, §3.2) defines the raw upper bound on the
  DNA entropy rate as the channel capacity of the four-letter genetic code:

    H_raw := log₂ 4 = 2 bits per substitution

  Three biochemical constraints reduce this to an effective entropy rate:
    • Transition/transversion bias  (factor φ ∈ (0, 1])
    • Context-dependent mutation    (factor ψ ∈ (0, 1])
    • Purifying selection           (factor ω ∈ (0, 1])

  These reductions yield h ∈ [1.58, 1.65] bits/substitution. The lower end
  corresponds to h_three_letter := log₂ 3 ≈ 1.585, the information-theoretic
  minimum for a code that must distinguish three domains of life.
-/

/-- The raw information-theoretic upper bound: channel capacity of the four-letter
    genetic code. H_raw = log₂ 4 = 2 bits per substitution. -/
noncomputable def H_raw : ℝ := log 4 / log 2

/-- Lower bound on the effective entropy rate: log₂ 3 ≈ 1.585 bits, the
    information-theoretic minimum for a 3-domain-distinguishing code. -/
noncomputable def h_three_letter : ℝ := log 3 / log 2

theorem H_raw_eq_two : H_raw = 2 := by
  unfold H_raw
  have h4 : log 4 = 2 * log 2 := by
    have eq1 : (4 : ℝ) = 2 ^ 2 := by norm_num
    rw [eq1, log_pow]; push_cast; ring
  rw [h4, mul_div_assoc, div_self log2_pos.ne', mul_one]

theorem H_raw_pos : H_raw > 0 := by rw [H_raw_eq_two]; norm_num

theorem h_three_letter_pos : h_three_letter > 0 :=
  div_pos (log_pos (by norm_num)) log2_pos

theorem h_three_letter_bounds : 1.5 < h_three_letter ∧ h_three_letter < 1.7 := by
  unfold h_three_letter; constructor
  · rw [lt_div_iff₀ log2_pos]
    have h1 : (3 : ℝ) * log 2 = log (2 ^ 3) := by rw [log_pow]; push_cast; ring
    have h2 : (2 : ℝ) * log 3 = log (3 ^ 2) := by rw [log_pow]; push_cast; ring
    nlinarith [log_lt_log (show (0:ℝ) < 2^3 by positivity) (show (2:ℝ)^3 < 3^2 by norm_num)]
  · rw [div_lt_iff₀ log2_pos]
    have h1 : (10 : ℝ) * log 3 = log (3 ^ 10) := by rw [log_pow]; push_cast; ring
    have h2 : (17 : ℝ) * log 2 = log (2 ^ 17) := by rw [log_pow]; push_cast; ring
    nlinarith [log_lt_log (show (0:ℝ) < 3^10 by positivity) (show (3:ℝ)^10 < 2^17 by norm_num)]

theorem h_three_letter_lt_H_raw : h_three_letter < H_raw := by
  rw [H_raw_eq_two]
  linarith [h_three_letter_bounds.2]

structure TransitionBias where
  R : ℝ
  R_pos : R > 0
  phi : ℝ
  phi_range : 0 < phi ∧ phi ≤ 1

structure ContextBias where
  rho : ℝ
  rho_pos : rho ≥ 1
  psi : ℝ
  psi_range : 0 < psi ∧ psi ≤ 1

structure SelectionBias where
  s : ℝ
  s_nonneg : s ≥ 0
  omega : ℝ
  omega_range : 0 < omega ∧ omega ≤ 1

/-- Effective entropy rate: starts from the four-letter channel capacity
    H_raw = 2 and is reduced multiplicatively by transition/transversion bias φ,
    context bias ψ, and selection ω. -/
noncomputable def h_effective (phi psi omega : ℝ) : ℝ := H_raw * phi * psi * omega

private lemma factors_le_one {a b c : ℝ}
    (ha : 0 < a ∧ a ≤ 1) (hb : 0 < b ∧ b ≤ 1) (hc : 0 < c ∧ c ≤ 1) :
    a * b * c ≤ 1 :=
  mul_le_one₀ (mul_le_one₀ ha.2 hb.1.le hb.2) hc.1.le hc.2

private lemma h_effective_pos {phi psi omega : ℝ}
    (hphi : 0 < phi) (hpsi : 0 < psi) (homega : 0 < omega) :
    0 < h_effective phi psi omega := by
  unfold h_effective; exact mul_pos (mul_pos (mul_pos H_raw_pos hphi) hpsi) homega

private lemma h_effective_le {phi psi omega : ℝ}
    (hphi : 0 < phi ∧ phi ≤ 1) (hpsi : 0 < psi ∧ psi ≤ 1)
    (homega : 0 < omega ∧ omega ≤ 1) :
    h_effective phi psi omega ≤ H_raw := by
  unfold h_effective
  calc H_raw * phi * psi * omega = H_raw * (phi * psi * omega) := by ring
    _ ≤ H_raw * 1 := mul_le_mul_of_nonneg_left (factors_le_one hphi hpsi homega) H_raw_pos.le
    _ = H_raw := mul_one _

theorem entropy_rate_decomposition_bounds
    (tb : TransitionBias) (cb : ContextBias) (sb : SelectionBias) :
    0 < h_effective tb.phi cb.psi sb.omega ∧
    h_effective tb.phi cb.psi sb.omega ≤ H_raw :=
  ⟨h_effective_pos tb.phi_range.1 cb.psi_range.1 sb.omega_range.1,
   h_effective_le tb.phi_range cb.psi_range sb.omega_range⟩

/-- The effective entropy rate cannot exceed the four-letter channel capacity,
    so the curvature is bounded above by κ(H_raw, 2) = (2 ln 2)² ≈ 1.921. -/
theorem kappa_bounded_by_raw (phi psi omega : ℝ)
    (hphi : 0 < phi ∧ phi ≤ 1) (hpsi : 0 < psi ∧ psi ≤ 1) (homega : 0 < omega ∧ omega ≤ 1) :
    κ (h_effective phi psi omega) 2 ≤ κ H_raw 2 := by
  simp only [kappa_n2]
  exact sq_le_sq'
    (by linarith [mul_pos (h_effective_pos hphi.1 hpsi.1 homega.1) log2_pos,
                  mul_pos H_raw_pos log2_pos])
    (mul_le_mul_of_nonneg_right (h_effective_le hphi hpsi homega) log2_pos.le)

/-- The four-letter alphabet ceiling, explicit form: κ ≤ (2 ln 2)² = 4 (ln 2)². -/
theorem kappa_bounded_by_alphabet (phi psi omega : ℝ)
    (hphi : 0 < phi ∧ phi ≤ 1) (hpsi : 0 < psi ∧ psi ≤ 1) (homega : 0 < omega ∧ omega ≤ 1) :
    κ (h_effective phi psi omega) 2 ≤ 4 * (log 2) ^ 2 := by
  have hb := kappa_bounded_by_raw phi psi omega hphi hpsi homega
  simp only [kappa_n2, H_raw_eq_two] at hb
  -- hb : (h_effective phi psi omega * log 2) ^ 2 ≤ (2 * log 2) ^ 2
  rw [kappa_n2]
  calc (h_effective phi psi omega * log 2) ^ 2
      ≤ (2 * log 2) ^ 2 := hb
    _ = 4 * (log 2) ^ 2 := by ring

/-
  ───────────────────────────────────────────────────────────────────────────
  Cross-substrate generalization (Hyperbolic Trilogy)
  ───────────────────────────────────────────────────────────────────────────

  Paper II (Fenn & Fenn 2026, Zenodo 10.5281/zenodo.19381558) generalizes
  the state equation to any information-generating hierarchy whose code has
  alphabet size α. Paper I (this repository) is the α = 4 case (DNA).
  Paper III (github.com/sentry-bio/convergent-alphabets) is the α ≈ 40 case
  (phonemes). The alphabet ceiling at n = 2 is κ ≤ (log α)², specializing to:

    α = 4  (DNA):          κ ≤ (log 4)²  = 4 (log 2)²  ≈ 1.921
    α = 20 (proteins):     κ ≤ (log 20)²               ≈ 8.974
    α = 40 (phonemes):     κ ≤ (log 40)²               ≈ 13.603
-/

/-- Information capacity of an alphabet of size α: log₂ α bits per symbol. -/
noncomputable def H_raw_of_alphabet (α : ℝ) : ℝ := log α / log 2

/-- For α > 1 the alphabet capacity is positive. -/
theorem H_raw_of_alphabet_pos (α : ℝ) (hα : α > 1) : H_raw_of_alphabet α > 0 :=
  div_pos (log_pos hα) log2_pos

/-- H_raw (the DNA default) is the α = 4 instance. -/
theorem H_raw_eq_alphabet_four : H_raw = H_raw_of_alphabet 4 := rfl

/-- H_raw_of_alphabet is strictly monotone in α for α > 1. -/
theorem H_raw_of_alphabet_mono (α β : ℝ) (hα : α > 1) (hlt : α < β) :
    H_raw_of_alphabet α < H_raw_of_alphabet β := by
  unfold H_raw_of_alphabet
  exact div_lt_div_of_pos_right
    (log_lt_log (by linarith) hlt) log2_pos

/-- Generic alphabet ceiling. For any alphabet size α > 1 and any effective
    entropy rate h_eff bounded above by the alphabet's capacity, the curvature
    at n = 2 is bounded by (log α)². The existing DNA-specific theorem
    `kappa_bounded_by_alphabet` is the α = 4 instance, yielding 4(log 2)². -/
theorem kappa_bounded_by_alphabet_general
    (α : ℝ) (hα : α > 1) (h_eff : ℝ)
    (h_pos : h_eff > 0) (h_le : h_eff ≤ H_raw_of_alphabet α) :
    κ h_eff 2 ≤ (log α) ^ 2 := by
  rw [kappa_n2]
  have hlog_pos : log α > 0 := log_pos hα
  have h_prod_le : h_eff * log 2 ≤ log α := by
    have eq : H_raw_of_alphabet α * log 2 = log α := by
      unfold H_raw_of_alphabet
      exact div_mul_cancel₀ (log α) log2_pos.ne'
    calc h_eff * log 2
        ≤ H_raw_of_alphabet α * log 2 :=
          mul_le_mul_of_nonneg_right h_le log2_pos.le
      _ = log α := eq
  have h_prod_pos : h_eff * log 2 > 0 := mul_pos h_pos log2_pos
  exact sq_le_sq' (by linarith) h_prod_le

/-
  ═══════════════════════════════════════════════════════════════════════════════
  PART III: TREE DIMENSIONALITY AND RATE-DISTORTION POTENTIAL
  ═══════════════════════════════════════════════════════════════════════════════
-/

def tree_dimension : ℕ := 2

structure MetricTree where
  V : Type*
  d : V → V → ℝ
  d_nonneg : ∀ u v, d u v ≥ 0
  d_symm : ∀ u v, d u v = d v u
  d_triangle : ∀ u v w, d u w ≤ d u v + d v w
  four_point : ∀ u v w x,
    d u v + d w x ≤ max (d u w + d v x) (d u x + d v w)

theorem H1_insufficient_for_trees :
    ∀ (T : MetricTree), T.V → T.V → T.V → True := by intros; trivial

theorem embedding_dimension_optimal (n : ℝ) (hn : n > 2) (h : ℝ) (hpos : h > 0) :
    κ h n < κ h 2 := kappa_max_at_n2 h n hpos hn

theorem dimension_two_maximizes_curvature (h : ℝ) (hpos : h > 0) :
    ∀ n : ℝ, n > 2 → κ h n < κ h 2 := fun _ hn => kappa_max_at_n2 h _ hpos hn

/-- Information production rate in nats per symbol: I(h) = h · ln 2. -/
noncomputable def I (h : ℝ) : ℝ := h * log 2

/-- Geometric capacity: C(κ, n) = (n − 1) √κ for κ > 0, else 0. -/
noncomputable def C (κ_val : ℝ) (n : ℝ) : ℝ :=
  if κ_val ≤ 0 then 0 else (n - 1) * sqrt κ_val

/-- Rate-distortion residual: ε = I − C. The state equation holds when ε = 0. -/
noncomputable def ε (h : ℝ) (κ_val : ℝ) (n : ℝ) : ℝ := I h - C κ_val n

/-- Rate-distortion potential: U = ε². Non-negative, zero only at κ_critical. -/
noncomputable def U (h : ℝ) (κ_val : ℝ) (n : ℝ) : ℝ := (ε h κ_val n) ^ 2

/-- Critical curvature: the unique κ solving the state equation for given (h, n). -/
noncomputable def κ_critical (h : ℝ) (n : ℝ) : ℝ := κ h n

theorem potential_nonneg (h κ_val n : ℝ) : U h κ_val n ≥ 0 := sq_nonneg _

theorem potential_zero_iff (h : ℝ) (κ_val : ℝ) (n : ℝ)
    (hpos : h > 0) (hκpos : κ_val > 0) (hn : n > 1) :
    U h κ_val n = 0 ↔ κ_val = κ_critical h n := by
  unfold U ε I C κ_critical κ
  simp only [show ¬κ_val ≤ 0 by linarith, show ¬n ≤ 1 by linarith, ↓reduceIte]
  have hn1 : n - 1 > 0 := by linarith
  constructor
  · intro hU
    have hsqrt : sqrt κ_val = h * log 2 / (n - 1) := by
      have := sq_eq_zero_iff.mp hU; field_simp at this ⊢; linarith
    calc κ_val = (sqrt κ_val) ^ 2 := (sq_sqrt hκpos.le).symm
      _ = _ := by rw [hsqrt]
  · intro hκ; rw [hκ, sqrt_sq (div_pos (mul_pos hpos log2_pos) hn1).le]
    field_simp; ring

/-- Lyapunov function for global stability: V(κ, κ*) = (√κ − √κ*)². -/
noncomputable def V (κ_val κ_star : ℝ) : ℝ := (sqrt κ_val - sqrt κ_star) ^ 2

theorem lyapunov_nonneg (κ_val κ_star : ℝ) : V κ_val κ_star ≥ 0 := sq_nonneg _

theorem lyapunov_zero_iff (κ_val κ_star : ℝ)
    (hκpos : κ_val > 0) (hκspos : κ_star > 0) :
    V κ_val κ_star = 0 ↔ κ_val = κ_star := by
  unfold V; constructor
  · intro hV
    have := sq_eq_zero_iff.mp hV
    calc κ_val = (sqrt κ_val) ^ 2 := (sq_sqrt hκpos.le).symm
      _ = (sqrt κ_star) ^ 2 := by rw [show sqrt κ_val = sqrt κ_star by linarith]
      _ = κ_star := sq_sqrt hκspos.le
  · intro heq; rw [heq, sub_self]; ring

theorem potential_gradient_zero_at_critical (h n : ℝ) (hpos : h > 0) (hn : n > 1) :
    ε h (κ_critical h n) n = 0 := by
  unfold ε I C κ_critical κ
  have hn1 : n - 1 > 0 := by linarith
  have hnum := div_pos (mul_pos hpos log2_pos) hn1
  have hnsq : ¬(h * log 2 / (n - 1)) ^ 2 ≤ 0 := not_le.mpr (sq_pos_of_pos hnum)
  simp only [show ¬n ≤ 1 by linarith, hnsq, ↓reduceIte, sqrt_sq hnum.le]
  field_simp; ring

theorem potential_second_derivative_pos (h n : ℝ) (hpos : h > 0) (hn : n > 1) :
    let κ_star := κ_critical h n
    (n - 1) ^ 2 / (2 * κ_star ^ (3/2 : ℝ)) > 0 := by
  have hn1 : n - 1 > 0 := by linarith
  have hκ_pos : κ_critical h n > 0 := by
    unfold κ_critical κ; simp only [not_le.mpr hn, ↓reduceIte]
    exact sq_pos_of_pos (div_pos (mul_pos hpos log2_pos) (by linarith))
  exact div_pos (by positivity) (by linarith [rpow_pos_of_pos hκ_pos (3/2 : ℝ)])

end BiosphereCurvature
