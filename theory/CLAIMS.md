# Claim & Artifact Registry

The whole-program map — layers, the two manuscripts, and their seams — is
[`PROGRAM.md`](PROGRAM.md). This file is the ledger it refers to.

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
| L1.1 | \(\beta\le c\,h_{\mathrm{pack}}\) for faithful finite-rate representations | THEOREM (Lean) | `theory/lean/ActiveGeometry/Packing.lean` (`faithful_representation_addressable`) |
| L1.2 | Finiteness holds in every proper metric host | THEOREM (Lean) | `theory/lean/ActiveGeometry/Packing.lean` (`hasFinitePacking_of_properSpace`) |
| L1.3 | Block address capacity \(=\) metric packing number (finite radius) | IDENTITY (Lean) | `theory/lean/ActiveGeometry/Packing.lean` (`exists_optimal_blockCode`, `card_le_packingCount`) |
| L1.4 | Asymptotic block identity \(C_{\mathrm{block}}=c\,h_{\mathrm{pack}}\) | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §4 (Theorem 4.2) |
| L1.5 | Polynomial-growth exclusion (Corollary 4.3) | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §4 |
| L1.6 | Constrained-capacity ladder: relational ≤ causal ≤ persistent ≤ block | THEOREM (ordering) / OPEN (rung values) | `theory/ADDRESSABILITY_KERNEL.md` §3 |
| L1.7 | Slack decomposition: block slack = relational tax + utilization slack | THEOREM (definitional) | `theory/MATHEMATICAL_SPINE.md` §5 |
| L1.8 | Weighted relational capacity of \(\mathbb H_\kappa^n\) is \(c(n-1)\sqrt\kappa\) | THEOREM (paper; cited lower bound) | `theory/RELATIONAL_CAPACITY_THEOREM.md`; `theory/MATHEMATICAL_SPINE.md` §4 (Theorem 4.4) |
| L1.9 | Exact unit-edge formulation of former Conjecture 4.4 | REFUTED for \(c<\varepsilon\); open for \(c\ge\varepsilon\) | `theory/RELATIONAL_CAPACITY_THEOREM.md` §1, §5 |

## Bridge — four-point classifier

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| B.1 | Four-point condition ⇔ tree metric; \(\delta=0\) ⇔ 0-hyperbolic | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 (Theorem 6.1) |
| B.2 | \(\delta\) classifies tree-ness, not curvature magnitude | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 |

## Layer IIa — curvature realization (host class)

*Which geometry hosts the data. The better-supported biological claim. Its
premise — isotropy — is asserted, not measured; the decisive test is E9.*

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L2a.1 | Isotropic hyperbolic host has \(h_{\mathrm{vol}}=(n-1)\sqrt\kappa\) | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §7 |
| L2a.2 | Curvature floor \(\kappa\ge(\beta/(c(n-1)))^2\) (realization + bound, no saturation) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`isotropic_curvature_at_least_floor`) |
| L2a.3 | Curvature genericity: hyperbolic is the generic homogeneous realization of the relational-exponential class | OPEN | `theory/MATHEMATICAL_SPINE.md` §7 (Conjecture 7.1) |
| L2a.4 | \(n=2\) is an embeddability floor for branching trees, not a fitted constant | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §6, §7 |
| L2a.5 | Isotropy is an asserted premise (`--assume-isotropic-hyperbolic`), not a measurement | HONESTY ITEM | `tools/addressability_meter.py`; test E9 |

## Layer IIb — saturation

*Whether a process fills its budget, yielding the state-equation equality. The
harder, less-supported claim; independent tests currently fail their kill lines.*

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L2b.1 | Saturation condition \(\eta=1\) (coordinate-free) | DEFINITION | `theory/MATHEMATICAL_SPINE.md` §5 |
| L2b.2 | State equation \(\bar\kappa^*=(h_{\mathrm{eff}}\ln2/(n-1))^2\) under saturation + isotropy | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`normalized_state_equation`) |
| L2b.3 | \(\bar\kappa=c^2\kappa\) is the unit-invariant curvature; raw formula needs gauge \(c=1\) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`normalized_curvature_scale_invariant`, `process_time_gauge`) |
| L2b.4 | Some hierarchy actually saturates | EMPIRICAL, open (independent tests currently fail) | `theory/DECISIVE_EXPERIMENTS.md` E3, E4 |
| L2b.5 | Near-capacity hyperbolic codes concentrate retained histories at \(d(o,f(v))\approx c\tau(v)\); deficit tradeoff \(e^{-(h\delta c-\eta)R}\) | THEOREM (finite count Lean; asymptotic paper) | `theory/RELATIONAL_CAPACITY_THEOREM.md` §5; `theory/lean/ActiveGeometry/Packing.lean` (`subball_fraction_le_packing_fraction`) |

## Instruments (status of tools, not claims about nature)

| # | Instrument | Status | Backing artifact |
|---|---|---|---|
| I.1 | Radial-rate meter (M2) | CERTIFIED on synthetic ground truth | `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.2 | Tree-defect meter (M4) | CERTIFIED on synthetic ground truth | `tools/addressability_meter.py`; `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.3 | Packing-entropy magnitude (M3) | NOT CERTIFIED as a magnitude estimator | `theory/DECISIVE_EXPERIMENTS.md` E1 |
| I.4 | Growth-class gate (exponential vs polynomial) | CERTIFIED (13/13); makes Corollary 4.3 enforceable | `theory/DECISIVE_EXPERIMENTS.md` E1 |
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
| E1 meter certification | I | run (M2/M4 pass, M3 magnitude fails, gate passes) | instruments legal |
| E8 boundary mapping | I | designed | premises load-bearing |
| E9 matched-capacity Euclidean vs hyperbolic (*E-alpha*) | IIa | **designed, unrun — highest value** | is hyperbolic forced |
| E2 equal-edge numerical achievability | IIa refinement | run (endpoint obstruction) | stronger synchronization subclass; Theorem 4.4 already settles host capacity |
| E7 reticulation intervention | IIa | designed | \(\delta\perp h_{\mathrm{pack}}\) in vivo |
| E3 barcoded lineages | IIb | designed | \(\eta\le1\) with a given tree |
| E4 mutation-rate intervention | IIb | designed | saturation as a response |
| E5 trained hierarchy | IIb | designed (small-scale runnable) | saturation cross-domain |
| E6 radius concentration | IIa+IIb | designed | necessary shell consequence of near-capacity hyperbolic coding |

Balance: IIa (host class) is comparatively well-supported by existing evidence
yet under-tested by intervention — E9, the decisive one, is unrun. IIb
(saturation) is heavily tested and currently fails its independent kill lines.
Closing E9 is the program's highest-leverage empirical step.

## What is explicitly not claimed

- No absolute curvature \(\kappa\) is a certified measurement (M3 magnitude is
  uncertified; embedding non-identifiability is documented).
- No zero-free-parameter cross-domain law is established; the state equation is
  a Layer IIb ideal awaiting saturation evidence and the remaining
  curvature-genericity conjecture.
- The hyperbolic host class (IIa) is not yet *forced* by evidence: isotropy is
  asserted, and E9 has not been run.
- Lean certifies algebra from definitions; it certifies neither that biology
  saturates nor either open conjecture.
