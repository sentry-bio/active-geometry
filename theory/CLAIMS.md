# Claim & Artifact Registry

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
- **OPEN** — a stated conjecture, not proved.
- **CONVENTION** — a gauge/coordinate choice, true by stipulation.
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
| L1.8 | Relational capacity of \(\mathbb H_\kappa^n\) is \(c(n-1)\sqrt\kappa\) (subcritical achievability) | OPEN | `theory/MATHEMATICAL_SPINE.md` §4 (Conjecture 4.4) |

## Bridge — four-point classifier

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| B.1 | Four-point condition ⇔ tree metric; \(\delta=0\) ⇔ 0-hyperbolic | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 (Theorem 6.1) |
| B.2 | \(\delta\) classifies tree-ness, not curvature magnitude | THEOREM (paper) | `theory/MATHEMATICAL_SPINE.md` §6 |

## Layer II — curvature realization

| # | Claim | Status | Backing artifact |
|---|---|---|---|
| L2.1 | Isotropic hyperbolic host has \(h_{\mathrm{vol}}=(n-1)\sqrt\kappa\) | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §7 |
| L2.2 | Curvature floor \(\kappa\ge(\beta/(c(n-1)))^2\) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`isotropic_curvature_at_least_floor`) |
| L2.3 | State equation \(\bar\kappa^*=(h_{\mathrm{eff}}\ln2/(n-1))^2\) under saturation + isotropy | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`normalized_state_equation`) |
| L2.4 | \(\bar\kappa=c^2\kappa\) is the unit-invariant curvature; raw formula needs gauge \(c=1\) | THEOREM (Lean algebra) | `theory/lean/ActiveGeometry/Addressability.lean` (`normalized_curvature_scale_invariant`, `process_time_gauge`) |
| L2.5 | Curvature genericity: hyperbolic is the generic homogeneous realization of the relational-exponential class | OPEN | `theory/MATHEMATICAL_SPINE.md` §7 (Conjecture 7.1) |
| L2.6 | \(n=2\) is an embeddability floor for branching trees, not a fitted constant | THEOREM (cited) | `theory/MATHEMATICAL_SPINE.md` §6, §7 |

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
| E.3 | Radius is accumulated information (the \(c\)-axis) | EMPIRICAL, open | `theory/DECISIVE_EXPERIMENTS.md` E6 |

## What is explicitly not claimed

- No absolute curvature \(\kappa\) is a certified measurement (M3 magnitude is
  uncertified; embedding non-identifiability is documented).
- No zero-free-parameter cross-domain law is established; the state equation is
  a Layer II ideal awaiting the two open conjectures and the empirical program.
- Lean certifies algebra from definitions; it certifies neither that biology
  saturates nor either open conjecture.
