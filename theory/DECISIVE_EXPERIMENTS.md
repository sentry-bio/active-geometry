# Decisive Experiments for the Addressability Limit

## What an experiment can and cannot establish

The theorems are not at stake. The packing converse and the block identity are
counting facts, machine-checked, and no measurement can strengthen or damage
them. Experiments decide the two questions mathematics cannot:

1. **Applicability.** Do natural systems instantiate the premises — retained
   distinguishable histories, fixed operational resolution, finite radial
   rate?
2. **Saturation.** Where the premises hold, does the system fill its budget
   (\(\eta\to 1\)), and is it *driven* there?

Decisiveness is ranked: **interventions** (turn a knob the theory names,
predict the response function) beat **pre-registered predictions** (state the
number before measuring) beat **calibrated observations** (fit after the
fact). Each experiment below names its layer, its independence firewall, its
predicted outcome, and its kill line. A protocol without a kill line is not an
experiment.

Every experiment reports the full vector of
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md) §11:

\[
(\beta,\ c,\ h_{\mathrm{pack}},\ \eta,\ \delta,\ n,\ \bar\kappa)
\]

with per-quantity estimator provenance. No quantity may appear on both sides
of a comparison through a shared input. Violations reclassify the result as
CIRCULAR regardless of numerical agreement.

The core set is five experiments. Two more (E6, E7) close remaining
identifications. Together they cover the instrument, the open theorem, the
independence firewall, saturation under intervention, and the cross-domain
bet.

---

## E1 — Meter certification on synthetic ground truth

**Layer.** Instrument validity. Precondition for everything below.

**System.** Generated processes on a lattice of known
\((\beta, c, h_{\mathrm{pack}}, \delta)\): \(b\)-ary trees in
\(\mathbb H^2_\kappa\) across \(\kappa\); reticulated variants at controlled
transfer fraction; Euclidean and product-space nulls; non-stationary rate
schedules.

**Independent estimators.** Each field estimator is scored only against the
axis it claims to read. Capacity estimators must be blind to reticulation;
defect estimators must be zero on pure trees of every curvature; Euclidean
nulls must yield \(\hat h_{\mathrm{pack}}\) at polynomial (zero exponential)
rate.

**Prediction.** Every field estimator recovers its own axis within stated
error and reads null on the others.

**Kill line.** Any estimator that misses synthetic truth beyond its stated
error, or responds to the orthogonal axis, is disqualified from field use.
Published numbers produced by a disqualified estimator are withdrawn, not
defended.

**Why decisive.** The prior failure mode of the program — an \(\mathbb H^2\)
stress estimator that missed synthetic targets by 34–109% while producing the
published biological curvatures — is exactly what E1 exists to catch before
the field. An uncertified meter cannot confirm anything.

**Status.** Partially done (`kernel_orthogonality`-class runs; ball-counting
estimator validated on \(b\)-ary trees). The full lattice, including
non-stationary schedules and the \(c\)-estimator, is open.

---

## E2 — Numerical achievability at fixed host geometry

**Layer.** The relational capacity conjecture
([`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md) Conjecture 4.4).

**System.** At fixed \((\kappa, n=2, \varepsilon, c)\), explicit relational
codes of \(b\)-ary trees at increasing depth \(R\): Sarkar-type constructions
and learned embeddings. Record realized rate \(\hat\beta(R)\) and the smallest
distortion \((D,K)\) at which the code remains admissible.

**Prediction if the conjecture is true.** \(\hat\beta \to c(n-1)\sqrt\kappa\)
with \(D, K\) bounded in \(R\).

**Prediction if false.** Holding \(\hat\beta\) above some
\(\beta_0 < c(n-1)\sqrt\kappa\) forces \(D(R)\to\infty\); the plateau
\(\beta_0\) estimates the relational capacity and the gap becomes a measured
host invariant.

**Kill line.** None. Both outcomes are informative. This experiment cannot
prove the conjecture; it tells us which half to spend proof effort on.

**Why decisive.** It is the cheapest probe of the one open mathematical
problem. A bounded-distortion approach to the packing bound would license a
serious attempt at the coding theorem; a forced-distortion plateau would
reframe the program around a new host invariant.

---

## E3 — Ground-truth genealogy: barcoded lineages

**Layer.** Saturation, with the inference circularity removed entirely.

**System.** CRISPR lineage-recording systems (GESTALT-class barcode arrays and
expressed-barcode variants). The true genealogy is *written into* the cells
as an accumulating barcode, independently of any sequence-inferred tree.
Single-cell transcriptomes or other phenotypes supply the representation
metric. Organoid or embryo reconstructions give a known process clock.

**Independent estimators.**

- \(\beta\): barcode-edit rate per cell division (the recording array is the
  clock; it is not the phenotype).
- \(c\): radial growth of the phenotype embedding per division.
- \(h_{\mathrm{pack}}\): packing growth of the phenotype distance matrix.
- \(\delta\): quartet defect of the phenotype metric, scored against the
  *barcode* tree.

No quantity is inferred from the same object as any other.

**Prediction.** \(\eta \le 1\) at every depth (the bound). If description
length is under selection in the reconstruction, \(\eta\) is high and stable
across independent reconstructions. \(\delta\) of the phenotype metric
against the barcode tree is the HGT/homoplasy reading with a ground-truth
reference.

**Kill line.** A reconstruction with independently measured
\(\beta > c\,h_{\mathrm{pack}}\) at fixed resolution would challenge the
bound itself (and therefore the instrument, first). Systematic
\(\eta \ll 1\) in systems claimed to be under description-length pressure
kills saturation for this class without touching the bound.

**Why decisive.** It is the phylogenetic program with the tree given, not
inferred. Every prior circularity diagnostic in this program traces to
shared-input trees. E3 deletes that input.

---

## E4 — Mutation-rate intervention

**Layer.** Saturation as a *response*, not a correlation.

**System.** A serial-passage microbial or viral population with a controllable
mutation rate: mutator strains, mutagens at titrated dose, or a directed-
evolution setup with an adjustable polymerase error rate. Sample at matched
generation counts across mutation-rate conditions. The representation is a
phenotype or k-mer metric that does not take mutation rate as an input.

**Independent estimators.** \(\beta\) from the imposed (or independently
assayed) substitution rate per generation. \(c\) and \(h_{\mathrm{pack}}\)
from the representation metric only. \(\delta\) from the same metric, as a
negative control: mutation rate is not reticulation.

**Prediction.**

\[
\frac{\partial(c\,h_{\mathrm{pack}})}{\partial\beta}
\;\approx\;
1
\quad\text{in the saturated regime,}
\]

with \(\eta\) remaining high as \(\beta\) is moved. The bound
\(\eta\le 1\) is never violated. \(\delta\) is approximately invariant under
the mutation-rate knob (orthogonality).

**Kill line.** Capacity that does not track \(\beta\) (flat
\(c\,h_{\mathrm{pack}}\) under a several-fold change in mutation rate) kills
the claim that the system is driven toward saturation. \(\eta>1\) at
certified estimators kills the instrument or the bound. A \(\delta\) that
tracks mutation rate kills the orthogonality reading of the four-point
object.

**Why decisive.** This is the highest-ranked experiment in the set. The theory
names a knob (\(\beta\)) and a response (\(c\,h_{\mathrm{pack}}\)). Turning
the knob and watching the response is what distinguishes a law from a fit.
Observational \(\eta\approx 1\) across clades can always be a shared-input
artifact; a dose-response cannot.

**Near-term proxy.** Time-stamped serially sampled viruses (influenza,
SARS-CoV-2) give \(\beta\) and \(c\) in the same physical clock without
fossils. They are not an intervention, but they are the cleanest
observational version of the same identification and should be run first as
a protocol rehearsal.

---

## E5 — Description-length pressure in a trained hierarchy

**Layer.** The cross-domain saturation bet, with every variable under
experimental control.

**System.** A generative model trained on synthetic hierarchical data whose
alphabet entropy \(h\) is an experimental parameter. The learned
representation is the host. Training is explicit description-length pressure
(a rate-distortion, \(\beta\)-VAE, or hierarchical contrastive objective).
Measure \(\eta\) across epochs.

**Independent estimators.** \(\beta\) from the known generating process (not
from the network). \(c\) and \(h_{\mathrm{pack}}\) from the representation
geometry. The loss family is a third experimental axis.

**Prediction.** Under a lossless or near-lossless objective, \(\eta(t)\)
rises and plateaus near 1. Under a lossy objective, \(\eta\) plateaus below 1
at a value ordered by the rate-distortion parameter. Euclidean-constrained
architectures cannot host retained exponential novelty at finite radial rate
(Corollary 4.3): either \(\beta\to 0\), \(c\to\infty\), or collapse.

**Kill line.** A lossless objective whose \(\eta\) plateaus well below 1 after
capacity and finite-size effects are controlled kills the claim that
description-length pressure drives saturation. Successful Euclidean hosting
of retained exponential novelty at finite \(c\) kills the polynomial-
exclusion corollary in this operationalization.

**Why decisive.** It is the eka-silicon row of the substrate table, run as an
experiment rather than an observation. No GTDB, no shared tree, no biology.
If \(\eta\) does not respond to a knob that *is* description length, the
cross-domain saturation claim is not a law of information-generating
hierarchies; it is at best a fact about some evolved lineages.

---

## E6 — Radius identification (the \(c\)-axis)

**Layer.** The remaining identification in the CCS stack: whether radius is
accumulated information.

**System.** Lineages in which time, sequence divergence, and functional
complexity come apart. Bradytelic ("living fossil") clades, with deep
divergence from a sister group but low morphological and genomic change;
equivalently, bursts of functional innovation at modest elapsed time.

**Independent estimators.** Process time from fossils or a molecular clock
that is *not* the embedding radius. Sequence divergence from an alignment
that is *not* the embedding. Functional complexity from an independent
phenotype or gene-content measure. Embedded radius from a representation
that takes none of these as input.

**Prediction.** If radius is information, embedded \(r\) tracks functional
complexity when the three clocks disagree, and does not track elapsed time.
MDL already predicts this (the E10 remark); E6 is that prediction, tested.

**Kill line.** Embedded radius that tracks elapsed time rather than
complexity, in a sample where the two are known to disagree, kills the
identification of \(c\) with an information rate. The bound is untouched;
the CCS advisory axis remains advisory.

**Why decisive.** The kernel treats \(c\) as logically independent. CCS
currently certifies \(\theta\) and leaves \(r\) advisory (cross-instrument
Spearman \(0.46\)). Until E6, the instrument grounds *which* history, not
*how much room per step*.

---

## E7 — Reticulation intervention (\(\delta\perp h_{\mathrm{pack}}\))

**Layer.** The Buneman–Gromov orthogonality, in vivo.

**System.** A microbial population with a controllable horizontal-transfer
rate: conjugative plasmids at titrated donor density, or a phage-mediated
transduction gradient. The substitution process is held as constant as the
apparatus allows.

**Independent estimators.** Transfer rate from a marked mobile element
(counts of acquired markers, not from the phenotype metric).
\(h_{\mathrm{pack}}\) from the chromosome k-mer or phenotype metric.
\(\delta\) from the same metric.

**Prediction.** \(\delta\) tracks transfer rate. \(h_{\mathrm{pack}}\) is
approximately invariant. The four-point object splits on living cells the
way `kernel_orthogonality.py` splits on synthetic grids.

**Kill line.** A \(\delta\) that is blind to a several-fold change in transfer
rate, or an \(h_{\mathrm{pack}}\) that tracks transfer as if it were
branching, kills the operational orthogonality claim for this class.

**Why decisive.** It is E1's orthogonality test with the lattice replaced by
an organism. A theory that says volume and thinness are independent
coordinates of any distance matrix is required to survive a knob that moves
only one of them.

---

## Ranked reading

| Rank | Experiment | What it decides | Class |
|---|---|---|---|
| 1 | E4 mutation-rate intervention | Is saturation a response? | intervention |
| 2 | E5 trained hierarchy | Is saturation cross-domain? | intervention |
| 3 | E3 barcoded lineages | Is \(\eta\) real when the tree is given? | ground-truth observation |
| 4 | E7 reticulation intervention | Is \(\delta\perp h_{\mathrm{pack}}\) in vivo? | intervention |
| 5 | E6 radius identification | Is \(c\) information? | pre-registered prediction |
| 6 | E2 numerical achievability | Which half of Conjecture 4.4 to prove? | numerical |
| 7 | E1 meter certification | Are the instruments legal? | calibration |

E1 is last in the table and first in time. Nothing above it is interpretable
until the meter recovers synthetic truth.

## What would constitute firm theoretical ground

Not any one of these. The bound is already firm. Firm *applicability* is:

- E1 passed (legal instruments);
- E3 showing \(\eta\le 1\) with a given tree (premises instantiated, bound
  respected in a living hierarchy);
- E4 showing capacity tracking \(\beta\) (saturation as a response, not a
  correlation).

Firm *generality* is E5 in addition: the same response under a knob that is
description length and is not biology.

Firm *relational geometry* is a proof or a clean numerical settlement of
Conjecture 4.4 (E2 as reconnaissance).

A single failed intervention does not unwind the theorems. It reclassifies
saturation from a law-like regularity to a domain fact, which is a legitimate
and fully specified outcome of this protocol.
