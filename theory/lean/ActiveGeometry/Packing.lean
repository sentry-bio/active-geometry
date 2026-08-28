/-
  The addressability limit
  ========================

  This file is the geometric kernel. It proves the convergent-rate packing
  theorem that mathematicians should cite for the ordinary-limit case:

      convergent_rate_addressability_limit :
        FaithfulRepresentation → r(R) → ∞ →
          (three rates converge) → Addressable β c h_pack

  in every proper metric host, with packing-number finiteness discharged.

  A faithful representation supplies a source type of histories and, at
  every generative depth `R`:

  * a finite, nonempty census of source histories;
  * an address map from those histories into the metric host;
  * injectivity and fixed-resolution separation of their addresses;
  * containment of the addresses in a ball of radius `r(R)`.

  Mathlib's `Metric.packingNumber` is the maximal cardinality, in extended
  naturals, of a separated subset of a metric set. Hence the represented count
  is bounded by the packing number of its containing ball at every depth.
  Conversely, when that number is finite, an exact finite block code attaining
  it exists. Thus operational finite-block address capacity and metric packing
  capacity are identical in every proper metric host.

  If the radii tend to infinity and the following three *ordinary* limits
  exist (Lean `Tendsto`; the paper's limsup statement is strictly more
  general):

      log represented-count / R  → β,
      r(R) / R                    → c,
      log packing-count / ρ       → h_pack,

  then

      β ≤ c · h_pack.

  The conclusion is `Capacity.Addressable`. No curvature, isotropy,
  saturation, dimension, alphabet, or physical dynamics enters this proof.
  Retention (`histories_monotone`) is recorded on a separate structure; the
  bound does not use it. The block-achievability theorem does not claim that
  optimal codebooks across depths are nested, causal, or preserve a source
  hierarchy's relational metric.

  Not formalized here (or anywhere in this library): the paper limsup
  theorem; Theorem 4.4 (Skenderi / weighted relational capacity of `ℍⁿ_κ`);
  Theorem 7.1 (Heintze isotropy / axiom A3).
-/

import ActiveGeometry.Capacity
import Mathlib.Topology.MetricSpace.CoveringNumbers
import Mathlib.Topology.MetricSpace.ProperSpace
import Mathlib.Topology.Order.OrderClosed

namespace ActiveGeometry.Packing

open Filter
open scoped ENNReal NNReal Topology

variable {M : Type*} [MetricSpace M]

/-! ### Packing numbers -/

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

/-! ### Finite-block identity -/

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

/-- Finite-radius concentration converse. If `t` is any subcollection of a
    nonempty separated codebook `s` and `t` lies in a smaller ball, then the
    fraction of codewords in `t` is bounded by that smaller ball's exact
    packing count divided by the full codebook count.

    Applying this with `t` the radially deficient histories is the
    machine-checked counting step of the radial-concentration theorem. The
    exponential packing envelope and asymptotic rate comparison remain
    paper-level inputs. -/
theorem subball_fraction_le_packing_fraction
    (o : M) (ε : ℝ≥0) (ρ : ℝ) (s t : Finset M)
    (hρ : 0 ≤ ρ) (hfinite : HasFinitePacking o ε)
    (hs_nonempty : s.Nonempty) (hts : t ⊆ s)
    (hsep : Metric.IsSeparated ε (s : Set M))
    (hcontained : (t : Set M) ⊆ Metric.closedBall o ρ) :
    (t.card : ℝ) / (s.card : ℝ) ≤
      (packingCount o ε ρ : ℝ) / (s.card : ℝ) := by
  have htsep : Metric.IsSeparated ε (t : Set M) :=
    hsep.mono (by
      intro x hx
      exact_mod_cast hts (by exact_mod_cast hx))
  have hcard :=
    card_le_packingCount o ε ρ t hρ hfinite htsep hcontained
  have hcard_real :
      (t.card : ℝ) ≤ (packingCount o ε ρ : ℝ) := by
    exact_mod_cast hcard
  have hs_pos : 0 < (s.card : ℝ) := by
    exact_mod_cast Finset.card_pos.mpr hs_nonempty
  exact (div_le_div_iff_of_pos_right hs_pos).2 hcard_real

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

/-! ### Histories and faithful address maps -/

/-- A finite history census with a depth-dependent address map into the host.
    Faithfulness has two parts: no two histories in the census share an
    address, and distinct addresses remain separated at resolution `ε`.
    Containment records the radial budget.

    This structure does *not* require retention. The packing bound is a fact
    about each depth's census; a process that overwrites rather than accretes
    is still subject to the count, but its growth rate is then not a rate of
    retained history. -/
structure FaithfulRepresentation (H : Type*) (o : M) (ε : ℝ≥0) where
  histories : ℕ → Finset H
  address : ℕ → H → M
  radius : ℕ → ℝ
  resolution_pos : 0 < ε
  histories_nonempty : ∀ R, (histories R).Nonempty
  address_injective : ∀ R, Set.InjOn (address R) (histories R : Set H)
  separated :
    ∀ R, Metric.IsSeparated ε (address R '' (histories R : Set H))
  contained :
    ∀ R, address R '' (histories R : Set H) ⊆ Metric.closedBall o (radius R)

/-- Source-history growth observed at depth `R`. It is a retained-history
    rate only when the source censuses are nested, as in
    `RetainedRepresentation`. -/
noncomputable def historyRate
    {H : Type*} {o : M} {ε : ℝ≥0}
    (rep : FaithfulRepresentation H o ε) (R : ℕ) : ℝ :=
  Real.log ((rep.histories R).card : ℝ) / (R : ℝ)

/-- Radial distance used per generative step at depth `R`. -/
noncomputable def radialRate
    {H : Type*} {o : M} {ε : ℝ≥0}
    (rep : FaithfulRepresentation H o ε) (R : ℕ) : ℝ :=
  rep.radius R / (R : ℝ)

/-- Finite-depth addressability: faithfully addressed source histories cannot
    outnumber the exact packing count of their containing ball. -/
theorem history_card_le_packingCount
    {H : Type*} (o : M) (ε : ℝ≥0) (rep : FaithfulRepresentation H o ε)
    (hfinite : HasFinitePacking o ε) (R : ℕ) (hρ : 0 ≤ rep.radius R) :
    (rep.histories R).card ≤ packingCount o ε (rep.radius R) := by
  have henat :
      (((rep.histories R).card : ℕ) : ℕ∞) ≤
        Metric.packingNumber ε (Metric.closedBall o (rep.radius R)) := by
    calc
      (((rep.histories R).card : ℕ) : ℕ∞)
          = (rep.histories R : Set H).encard := by simp
      _ = (rep.address R '' (rep.histories R : Set H)).encard :=
        (rep.address_injective R).encard_image.symm
      _ ≤ Metric.packingNumber ε (Metric.closedBall o (rep.radius R)) :=
        (rep.separated R).encard_le_packingNumber (rep.contained R)
  have htop :
      Metric.packingNumber ε (Metric.closedBall o (rep.radius R)) ≠ ⊤ :=
    hfinite hρ
  simpa [packingCount] using ENat.toNat_le_toNat henat htop

/-- The finite-depth count inequality becomes a pointwise inequality between
    normalized logarithmic rates away from the irrelevant depth `R = 0`. -/
theorem historyRate_le_capacity_eventually
    {H : Type*} (o : M) (ε : ℝ≥0) (rep : FaithfulRepresentation H o ε)
    (hfinite : HasFinitePacking o ε)
    (hradius : Tendsto rep.radius atTop atTop) :
    ∀ᶠ R in atTop,
      historyRate rep R ≤
        radialRate rep R * packingRate o ε (rep.radius R) := by
  have hradius_pos : ∀ᶠ R in atTop, 0 < rep.radius R :=
    hradius (eventually_gt_atTop 0)
  filter_upwards [eventually_gt_atTop 0, hradius_pos] with R hR hr
  have hcount :=
    history_card_le_packingCount o ε rep hfinite R hr.le
  have hcard_pos : 0 < ((rep.histories R).card : ℝ) := by
    exact_mod_cast Finset.card_pos.mpr (rep.histories_nonempty R)
  have hcount_real :
      ((rep.histories R).card : ℝ) ≤
        (packingCount o ε (rep.radius R) : ℝ) := by
    exact_mod_cast hcount
  have hlog :
      Real.log ((rep.histories R).card : ℝ) ≤
        Real.log (packingCount o ε (rep.radius R) : ℝ) :=
    Real.log_le_log hcard_pos hcount_real
  have hR_real : 0 < (R : ℝ) := by exact_mod_cast hR
  unfold historyRate radialRate packingRate
  calc
    Real.log ((rep.histories R).card : ℝ) / (R : ℝ)
        ≤ Real.log (packingCount o ε (rep.radius R) : ℝ) / (R : ℝ) :=
      (div_le_div_iff_of_pos_right hR_real).2 hlog
    _ = (rep.radius R / (R : ℝ)) *
          (Real.log (packingCount o ε (rep.radius R) : ℝ) / rep.radius R) := by
      field_simp [hR_real.ne', hr.ne']

/-- Convergent-rate addressability limit, with packing finiteness as an
    explicit hypothesis. The three rate hypotheses are independent:
    source-history growth, radial calibration, and host packing growth
    are not defined from one another. Radius divergence is stated separately
    because it is an asymptotic sampling premise, not part of faithfulness.

    The paper's theorem uses limsup; this declaration uses ordinary finite
    limits (`Tendsto`). -/
theorem convergent_rate_addressability_limit_of_hasFinitePacking
    {H : Type*} (o : M) (ε : ℝ≥0) (rep : FaithfulRepresentation H o ε)
    (hfinite : HasFinitePacking o ε)
    (β c hpack : ℝ)
    (hradius : Tendsto rep.radius atTop atTop)
    (hdemand : Tendsto (historyRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 hpack)) :
    Capacity.Addressable β c hpack := by
  have hsampled :
      Tendsto (fun R ↦ packingRate o ε (rep.radius R)) atTop (𝓝 hpack) := by
    simpa [Function.comp_def] using hpacking.comp hradius
  apply le_of_tendsto_of_tendsto hdemand (hradial.mul hsampled)
  exact historyRate_le_capacity_eventually o ε rep hfinite hradius

/-- The convergent-rate corollary of the Addressability Limit in a proper
    metric host.

    If a faithful finite-resolution representation has radii tending to
    infinity, convergent source-history growth `β`, convergent radial rate `c`,
    and the host has convergent packing growth `h_pack`, then
    `β ≤ c · h_pack`. Finiteness of packing numbers is a theorem for this
    host class, not an extra assumption.

    Cite this declaration for the Lean-checked ordinary-limit case. The
    paper's limsup formulation is the Addressability Limit. Theorem 4.4
    (weighted relational capacity of `ℍⁿ_κ`) and Theorem 7.1 (Heintze / A3)
    are not this theorem and are not Lean. -/
theorem convergent_rate_addressability_limit
    {H : Type*} [ProperSpace M]
    (o : M) (ε : ℝ≥0) (rep : FaithfulRepresentation H o ε)
    (β c hpack : ℝ)
    (hradius : Tendsto rep.radius atTop atTop)
    (hdemand : Tendsto (historyRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 hpack)) :
    Capacity.Addressable β c hpack :=
  convergent_rate_addressability_limit_of_hasFinitePacking o ε rep
    (hasFinitePacking_of_properSpace o rep.resolution_pos.ne') β c hpack
    hradius hdemand hradial hpacking

/-- A zero-capacity host cannot carry positive source-history growth under the
    hypotheses of the packing theorem. -/
theorem no_positive_growth_at_zero_capacity
    {H : Type*} (o : M) (ε : ℝ≥0) (rep : FaithfulRepresentation H o ε)
    (hfinite : HasFinitePacking o ε)
    (β c : ℝ) (hβ : 0 < β)
    (hradius : Tendsto rep.radius atTop atTop)
    (hdemand : Tendsto (historyRate rep) atTop (𝓝 β))
    (hradial : Tendsto (radialRate rep) atTop (𝓝 c))
    (hpacking : Tendsto (packingRate o ε) atTop (𝓝 0)) :
    False := by
  have hbound :=
    convergent_rate_addressability_limit_of_hasFinitePacking
      o ε rep hfinite β c 0 hradius hdemand hradial hpacking
  unfold Capacity.Addressable at hbound
  nlinarith

/-! ### Retention in the source -/

/-- A faithful representation whose source census retains prior histories.
    Addresses may change with depth; retention constrains histories, not host
    coordinates. The packing bound itself needs only faithfulness at each
    depth. -/
structure RetainedRepresentation (H : Type*) (o : M) (ε : ℝ≥0)
    extends FaithfulRepresentation H o ε where
  histories_monotone : ∀ R, histories R ⊆ histories (R + 1)

/-- Retention makes source-history counts nondecreasing in depth. -/
theorem history_card_mono {H : Type*} {o : M} {ε : ℝ≥0}
    (rep : RetainedRepresentation H o ε) (R : ℕ) :
    (rep.histories R).card ≤ (rep.histories (R + 1)).card :=
  Finset.card_le_card (rep.histories_monotone R)

end ActiveGeometry.Packing
