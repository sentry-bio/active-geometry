# Claim & Artifact Registry

The minimal dependency graph is [`THROUGHLINE.md`](THROUGHLINE.md). The
whole-program administrative map — layers, manuscripts, and seams — is
[`PROGRAM.md`](PROGRAM.md). This file is the ledger they refer to.

A single machine-checkable ledger of the program's load-bearing claims, each
tagged by status and bound to the artifact that backs it. This exists so that
documentation cannot outrun code again: every backticked repository path in a
tracked Markdown file is verified to exist on disk by
[`tools/check_doc_artifacts.py`](../tools/check_doc_artifacts.py), and every
claim below names where it is discharged.

The checker enforces the theory program's directories (`theory`, `tools`,
`tests`) as a fatal gate; `--all` scans the whole repository as a non-fatal
advisory (it currently surfaces pre-existing data/figure references in the
biology docs, a separate cleanup). Wire the enforced call into CI.

Status vocabulary:

- **THEOREM** — proved; Lean-checked where marked, else a paper proof.
- **IDENTITY** — proved equality of two independently defined quantities.
- **DEFINITION** — fixes vocabulary; carries no empirical content.
- **OPEN** — a stated conjecture, not proved.
- **REFUTED** — false under its stated quantifiers; counterexample recorded.
- **CONVENTION** — a gauge/coordinate choice, true by stipulation.
- **HONESTY ITEM** — a premise or limitation that must accompany the claim.
- **INSTRUMENT** — the status of a measurement tool, not a claim about nature.
- **EMPIRICAL** — a claim about nature; decided only by experiment.

## Layer I — universal capacity theory (curvature-free)

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L1.1 | \(\beta\le c\,h_{\mathrm{pack}}\) for finite source censuses with injective, \(\varepsilon\)-separated address maps (diverging radii, convergent rates, proper host) | THEOREM (Lean convergent-rate corollary; limsup in paper) | `theory/lean/ActiveGeometry/Packing.lean` (`convergent_rate_addressability_limit`) |
| L1.2 | Finiteness holds in every proper metric host | THEOREM (Lean) | `theory/lean/ActiveGeometry/Packing.lean` (`hasFinitePacking_of_properSpace`) |
| L1.3 | Block address capacity \(=\) metric packing number (finite radius) | IDENTITY (Lean) | `theory/lean/ActiveGeometry/Packing.lean` (`exists_optimal_blockCode`, `card_le_packingCount`) |
| L1.4 | Asymptotic block identity \(C_{\mathrm{block}}=c\,h_{\mathrm{pack}}\) | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §4 (Theorem 4.2) |
| L1.5 | Polynomial-growth exclusion / trichotomy (Corollary 4.3): exponential room, exponential addresses, or forgetting | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §4 |
| L1.6 | Constrained-capacity ladder: nested admissibility classes satisfy relational ≤ causal ≤ persistent ≤ block | DEFINITION (nested classes) / OPEN (rung values except the block ceiling) | `theory/ADDRESSABILITY_KERNEL.md` §3 |
| L1.7 | Slack decomposition: block slack = relational tax + utilization slack | THEOREM (definitional) | `theory/MATHEMATICAL_SPINE.md` §5 |

## Host-class closure of the relational rung

*Not curvature-free Layer I. This is the weighted relational identity in a real hyperbolic host; it closes the ladder's bottom rung in that host class.*

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L1.8 | Weighted relational capacity of \(\mathbb H_\kappa^n\) is \(c(n-1)\sqrt\kappa\) (supremum, not endpoint attainment) | THEOREM (paper; cited lower bound; host-specific) | `theory/RELATIONAL_CAPACITY_THEOREM.md`; `theory/MATHEMATICAL_SPINE.md` §4 (Theorem 4.4) |
| L1.9 | Exact unit-edge formulation of former Conjecture 4.4 | REFUTED for \(c<\varepsilon\); open for \(c\ge\varepsilon\) | `theory/RELATIONAL_CAPACITY_THEOREM.md` §1, §5 |

## Bridge — four-point classifier

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| B.1 | Four-point condition ⇔ tree metric; \(\delta=0\) ⇔ 0-hyperbolic | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 (Theorem 6.1) |
| B.2 | \(\delta\) classifies tree-ness, not curvature magnitude | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 |

## Layer IIa — curvature realization (host class)

*Which geometry hosts the data. The better-supported biological claim is
occupancy of exponential tree-like room. Isotropy is asserted, not
measured. E9 is a finite-sample shadow of Corollary 4.3, not a forcing
test.*

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L2a.1 | Isotropic hyperbolic host has \(h_{\mathrm{vol}}=(n-1)\sqrt\kappa\) | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §7 |
| L2a.2 | Curvature floor \(\kappa\ge(\beta/(c(n-1)))^2\) (space-form identification + bound, no saturation) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Capacity.lean` (`addressable_spaceForm_floor`) |
| L2a.3 | Inside the Heintze class, axiom A3 (full \(O(d)\) directional symmetry as host automorphisms) forces real \(\mathbb H^{d+1}\) | THEOREM (paper sketch; conditional on A3) | `theory/MATHEMATICAL_SPINE.md` §7 (Theorem 7.1) |
| L2a.4 | \(n=2\) is an embeddability floor for branching trees, not a fitted constant | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §6, §7 |
| L2a.5 | Isotropy is an asserted premise (`--assume-isotropic-hyperbolic`), not a measurement | HONESTY ITEM | `tools/addressability_meter.py` |
| L2a.6 | Axiom A3 holds of a real generator (no privileged directional structure) | HONESTY ITEM / EMPIRICAL, open | `theory/MATHEMATICAL_SPINE.md` §7; sector-wise \(h_{\mathrm{pack}}\) |

## Layer IIb — saturation

*Whether a process fills its budget, yielding the state-equation equality. The
harder, less-supported claim; independent tests currently fail their kill lines.*

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L2b.1 | Saturation condition \(\eta=1\) (coordinate-free) | DEFINITION | `theory/MATHEMATICAL_SPINE.md` §5 |
| L2b.2 | State equation \(\bar\kappa^*=(h_{\mathrm{eff}}\ln2/(n-1))^2\) under saturation + space-form identification | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/StateEquation.lean` (`normalized_state_equation`) |
| L2b.3 | \(\bar\kappa=c^2\kappa\) is the unit-invariant curvature; raw formula needs gauge \(c=1\) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Capacity.lean` (`normalized_curvature_scale_invariant`, `process_time_gauge`) |
| L2b.4 | Some hierarchy actually saturates | EMPIRICAL, open (independent tests currently fail) | `theory/DECISIVE_EXPERIMENTS.md` E3, E4 |
| L2b.5 | Near-capacity hyperbolic codes concentrate clock-relative radius; additive-deficit tradeoff \(e^{-(h\delta c-\Delta_{\rm cap})R}\), \(\Delta_{\rm cap}=ch(1-\eta)\) | THEOREM (finite count Lean; asymptotic paper with limsup qualification) | `theory/RELATIONAL_CAPACITY_THEOREM.md` §5; `theory/lean/ActiveGeometry/Packing.lean` (`subball_fraction_le_packing_fraction`) |

## Instruments (status of tools, not claims about nature)

| # | Instrument | Status | Backing artifact |
|---|---|---|---|
| I.1 | Radial-rate meter (M2) | CERTIFIED on synthetic ground truth | `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.2 | Tree-defect meter (M4) | CERTIFIED on synthetic ground truth | `tools/addressability_meter.py`; `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.3 | Packing-entropy magnitude (M3) | NOT CERTIFIED as a magnitude estimator | `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.4 | Growth-class gate (exponential vs polynomial) | CERTIFIED on full-span synthetics (13/13 in E1); refuses short radial windows | `tools/growth_class_gate.py`; `theory/MEASURABILITY.md`; `theory/lean/ActiveGeometry/Measurability.lean`; `tests/test_growth_class_gate.py` |
| I.5 | Independence firewall | CONVENTION only; runtime provenance check not yet implemented | `tools/addressability_meter.py` (`independence.verified: false`) |
| I.6 | Meter refuses to back-solve a missing axis from curvature | INSTRUMENT (enforced) | `tools/addressability_meter.py`; `tests/test_addressability_meter.py` |

## Empirical claims (decided only by experiment)

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| E.1 | Some hierarchy saturates its budget (\(\eta\to1\) / high utilization) | EMPIRICAL, open | `theory/DECISIVE_EXPERIMENTS.md` E3, E4 |
| E.2 | The law extends beyond biology | EMPIRICAL, open | `theory/DECISIVE_EXPERIMENTS.md` E5 |
| E.3 | Biological radius concentrates at accumulated process duration | EMPIRICAL, open; conditional on IIa+IIb | `theory/DECISIVE_EXPERIMENTS.md` E6 |

## Experiment allocation (which layer each test interrogates)

The biological claim is the host class (IIa); the protocol must not aim only at
saturation (IIb). This table is the allocation of record; see
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md) for the full designs.

| Experiment | Layer | Status | Decides |
|---|---|---|---|
| E1 meter certification | I | run (M2/M4 pass, M3 magnitude fails, gate passes on full-span synthetics) | instruments legal |
| E8 boundary mapping | I | designed | premises load-bearing |
| E9 matched-capacity Euclidean vs hyperbolic (*E-alpha*) | I shadow | designed, optional | is Corollary 4.3 visible at finite depth, matched packing (not "is hyperbolic forced") |
| E2 equal-edge numerical achievability | IIa refinement | run (endpoint obstruction) | stronger synchronization subclass; Theorem 4.4 already settles host capacity |
| E7 reticulation intervention | IIa | designed | \(\delta\perp h_{\mathrm{pack}}\) in vivo |
| E3 barcoded lineages | IIb | designed | \(\eta\le1\) with a given tree |
| E4 mutation-rate intervention | IIb | designed | saturation as a response |
| E5 trained hierarchy | IIb | rehearsal stopped: embedder lacks exponential fan-out ([`experiments/E5_DESIGN_FLAW.md`](../experiments/E5_DESIGN_FLAW.md)) | saturation cross-domain |
| E6 radius concentration | IIa+IIb | designed | necessary shell consequence of near-capacity hyperbolic coding |

Balance: IIa (occupancy) is comparatively well-supported by existing evidence
and is Paper II's figure. E9 is not the program's highest-leverage step; it
is an optional shadow of a theorem. IIb (saturation) is heavily tested and
currently fails its independent kill lines.

## What is explicitly not claimed

- No absolute curvature \(\kappa\) is a certified measurement (M3 magnitude is
  uncertified; embedding non-identifiability is documented).
- No zero-free-parameter cross-domain law is established; the state equation is
  a Layer IIb ideal awaiting saturation evidence
  ([`IIB_CONTRACT.md`](IIB_CONTRACT.md)). Former Conjecture 7.1 is
  Theorem 7.1, conditional on axiom A3; whether A3 holds of a real generator
  is empirical (L2a.6).
- The hyperbolic host class (IIa) is occupancy of exponential tree-like
  room, not a bake-off: isotropy is asserted; Corollary 4.3 already
  excludes polynomial hosts at finite rate; E9 does not force \(\mathbb H^2\).
- Lean certifies algebra from definitions; it certifies neither that biology
  saturates nor that axiom A3 holds, and it does not check Theorem 7.1.
