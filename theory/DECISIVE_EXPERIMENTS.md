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
fact).

This document is written to be executed. Every experiment states its inputs,
a numbered procedure, the estimator it calls (from the shared library below),
and a **decision rule** with explicit numeric thresholds — a predicted outcome
and a kill line. A protocol without a kill line is not an experiment.

The set is eight experiments: E1 (instrument), E2 (open theorem), E3–E4 and
E7–E8 (biology, mostly interventions), E5 (cross-domain), E6 (the \(c\)-axis).

---

## Shared conventions

**Pre-registration.** Before any field data is touched, write the predicted
number, the kill threshold, the random seeds, and the estimator version to a
frozen `preregistration.json`. Results computed by code that postdates that
file are calibrated observations, not predictions, and are labelled as such.

**The report vector.** Every run emits the
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md) §11 vector

\[
(\beta,\ c,\ h_{\mathrm{pack}},\ \eta,\ \delta,\ n,\ \bar\kappa)
\]

with a bootstrap confidence interval and a provenance tag per quantity naming
the object it was computed from (`clock`, `phenotype_metric`, `barcode_tree`,
`generator`, …).

**The independence firewall (the one hard rule).** No two quantities in a
comparison may share a provenance tag. \(\beta\) and \(h_{\mathrm{pack}}\)
compared through the same distance matrix is CIRCULAR and is reported as such
regardless of numeric agreement. Every experiment below is designed around
this rule.

*Enforcement status (corrected).* The rule is **not** yet mechanically
enforced. The meter (`tools/addressability_meter.py`) enforces weaker
disciplines: it refuses to silently promote the occupancy slope to host
entropy, refuses to compute \(\eta\) unless all three axes are supplied
independently on the command line, and always reports
`independence.verified: false` with a note that it cannot prove the axes were
estimated from different objects. A runtime provenance-tag check that refuses
\(\eta\) on tag collision is a **recommended addition, not present today**. Until
it exists, the firewall is a convention the operator must uphold, and the
meter's own disclaimer says so.

**Resolution.** Fix \(\varepsilon\) once per system and hold it across all
conditions. If distinguishability sharpens with depth, model
\(\varepsilon(R)\) explicitly (spine §1); do not let it float silently.

---

## Shared measurement procedures

These are the only estimators the experiments call. Each is certified once in
E1 and then reused. "Certified" means it passed E1's recovery and orthogonality
tests at its stated error.

- **M1 — \(\beta\) (retained-history rate).** From a *process clock that is not
  the representation*: barcode edits per division, assayed substitutions per
  generation, generating-process branch count, sound-change transitions. Fit
  \(\log N(R)\) vs depth \(R\); \(\beta\) is the slope. Provenance tag: `clock`.

- **M2 — \(c\) (radial rate).** Slope of representation radius \(r(R)\) from a
  fixed origin vs process depth \(R\). Origin is a declared anchor. Provenance
  tag: `embedding_radius`.

- **M3 — \(h_{\mathrm{pack}}\) (host capacity).** Ball-occupancy slope: for the
  representation distance matrix, count \(\varepsilon\)-separated points inside
  radius \(\rho\), regress \(\log P(\rho,\varepsilon)\) on \(\rho\). This is the
  meter's certified occupancy-slope mode, **not** a hyperbolic-stress fit.
  Provenance tag: `phenotype_metric` (or `representation_metric`).

- **M4 — \(\delta\) (tree defect).** Buneman four-point / Gromov defect over
  sampled quartets of the representation metric. Reported as the fraction of
  quartets with normalized slack above threshold. Provenance tag: same metric
  as M3 (so \(\delta\) and \(h_{\mathrm{pack}}\) are two readings of one matrix
  — that is allowed; they are never divided into each other).

- **M5 — \(\eta\) and \(\bar\kappa\).** \(\eta = \beta / (c\,h_{\mathrm{pack}})\),
  computed only when M1 and M3 carry different provenance tags. \(\bar\kappa\)
  and \(n\) only when an isotropic-hyperbolic fit is independently justified
  (spine §7); otherwise report \(h_{\mathrm{pack}}\) as the state variable and
  leave \(\bar\kappa\) blank.

The meter (`tools/addressability_meter.py`) already implements M3–M5 with the
refusal rule. M1 and M2 are supplied per system by the experiments below.

---

## E1 — Meter certification on synthetic ground truth

**Question.** Are the instruments legal? Precondition for everything else.

**Inputs.** A generator producing labelled distance matrices with known
\((\beta, c, h_{\mathrm{pack}}, \delta)\).

**Procedure.**
1. Generate a lattice: \(b\)-ary trees (\(b\in\{2,3,4,5\}\)) embedded in
   \(\mathbb H^2_\kappa\) across \(\kappa\in\{0.5,1,2,4\}\); reticulated
   variants at transfer fraction \(f\in\{0,0.01,0.05,0.25\}\); Euclidean and
   product-space nulls; two non-stationary rate schedules (rate step, rate
   ramp).
2. Subsample every matrix to a fixed leaf count, so a per-axis effect cannot
   be a leaf-count artifact.
3. Run M2, M3, M4 on every matrix. Run M1 on the generator.
4. Regress each estimator against every ground-truth axis, reporting both
   correlation **and absolute error** (correlation is scale-blind — a
   near-unity \(r\) can coexist with a large magnitude error).

**Decision rule.**
- *Pass:* an estimator recovers its own axis within a stated **absolute** error
  and reads null on the others.
- *Kill:* any estimator that misses its own axis beyond stated error, or moves
  with the orthogonal axis, is **disqualified from field use**, and any
  published number produced by it is withdrawn.

**Status — executed (2026-08-13).**
- **M2 (radial rate \(c\)): CERTIFIED.** Exact recovery, 0.00% error across all
  16 cells, \(r=1.0000\).
- **M4 (tree defect \(\delta\)): CERTIFIED.** Exactly zero on pure trees at
  every branching factor, tracks reticulation at \(r=0.946\), orthogonal to
  branching. This also reproduces the \(\delta\perp h_{\mathrm{pack}}\) split
  on synthetic ground truth.
- **M3 (packing entropy): NOT CERTIFIED as a single-point magnitude
  estimator.** On an Arm-A grid built so \(\eta=1\) exactly, measured \(\eta\)
  ran 0.87–1.19 with 11/16 cells above 1 — apparent bound violations from
  estimator noise alone. It fails the Euclidean null only slowly and breaks
  under non-stationarity (a branching-rate step gives ~25% error, driving
  \(\eta\) to 0.80 on data that sits exactly at 1).

**The growth-class gate (delivered and validated).** Exponential growth is
linear in \(\log P\) vs \(\rho\); polynomial growth is linear in \(\log P\) vs
\(\log\rho\). Fitting both and comparing adjusted \(R^2\) classified 13/13 hosts
correctly with no overlap (trees and \(\mathbb H^2\) exponential; 2D/3D grids
and Euclidean MSTs polynomial). **Certified M3 must carry this gate**: report a
packing entropy only when the exponential model wins, else \(h_{\mathrm{pack}}=0\)
by Corollary 4.3. This makes polynomial-growth exclusion — the coordinate-free
core of the theory — mechanically enforceable, and is the part of M3 that is
usable today. A **minimum-radial-shells precondition** is also required: M3
error tracks the number of fit points directly.

**Consequence for the protocol.** The *qualitative* axis (is the host
exponential?) is instrumented now. The *quantitative* axis (what is \(\eta\)?)
is not: an \(\eta\) within roughly \(\pm20\%\) of 1 is currently uninterpretable,
which is the entire decision range of E3 (\(\eta\ge0.8\)) and E5
(\(\eta\to1\) vs controls). E1 is correctly the gate on everything else, and it
does not yet fully pass.

---

## E2 — Numerical achievability at fixed host geometry

**Question.** Which half of Conjecture 4.4 deserves proof effort?

**Inputs.** Fixed \((\kappa, n=2, \varepsilon, c)\); a tree family of growing
depth \(R\).

**Procedure.**
1. For each depth \(R\), construct relational codes of a \(b\)-ary tree into
   \(\mathbb H^2_\kappa\). **Pin the construction**: children of a node must
   **subdivide their parent's angular sector** (a fixed cone half-angle is not
   admissible — it fails to partition angular space among children and yields
   an artifactual verdict). Run (a) the deterministic subdivided-sector layout
   and (b) a gradient-trained embedding minimizing distortion at fixed radius
   budget \(cR\).
2. Record realized growth \(\hat\beta(R)\) and the smallest \((D,K)\) for which
   conditions 1–4 of Conjecture 4.4 still hold.
3. Plot \(\hat\beta(R)\) against the bound \(c(n-1)\sqrt\kappa\), and
   \((D,K)\) against \(R\).

**Decision rule.**
- *Conjecture-true signature:* \(\hat\beta(R)\to c(n-1)\sqrt\kappa\) with
  \((D,K)\) bounded in \(R\) → attempt the coding theorem.
- *Conjecture-false signature:* maintaining \(\hat\beta\ge\beta_0\) forces
  \(D(R)\to\infty\); \(\beta_0\) estimates \(C_{\rm rel}\) and
  \(c(n-1)\sqrt\kappa-\beta_0\) is a new host invariant → redirect proof effort
  to the gap.
- *No kill line:* both outcomes are informative; this is reconnaissance for a
  proof, not a test of nature.

**Status — executed (2026-08-13), subdivided-sector construction.** At the
saturating rate \(c=\ln(b)/\sqrt\kappa\), conditions 1–3 hold exactly (radius =
budget to four decimals; separation exactly constant, \(d(\log\varepsilon)/dR
= -0.00000\), because arc length per node \(\sim\pi(e^\tau/b)^R\) is constant
exactly when \(\tau=\ln b\)). **Condition 4 fails linearly:** \(D(R)=0.558\,R
+ 0.034\), \(R^2=0.99996\) (a linear fit beats exponential). The binding pair
is angularly adjacent but tree-distant — neighbours on the circle whose
subtrees split at the root — which is exactly the §4.2 counterexample
mechanism made concrete. Extra radius does not rescue it: the distortion-optimal
rate sits stably at \(1.20\times\) saturating across depths, and minimum \(D\)
still grows (slope only drops \(0.56\to0.25\)).

This is the **conjecture-false signature at the saturating rate**, from two
constructions — reconnaissance, not proof. It redirects effort toward proving a
**strict gap** \(C_{\rm rel}<c(n-1)\sqrt\kappa\) with a lower bound on
distortion growth, i.e. a genuine price of genealogy, rather than toward proving
equality. A prior fixed-cone implementation reached the opposite (also
"false", but by an artifactual \(\varepsilon\)-decay mechanism); that two
implementers of the same prose diverged is why step 1 now pins the
construction.

---

## E3 — Ground-truth genealogy: barcoded lineages

**Question.** Is \(\eta\) real when the tree is given rather than inferred?

**Inputs.** A CRISPR lineage-recording dataset (GESTALT / scGESTALT /
expressed-barcode; zebrafish, mouse, or organoid) with (a) an accumulating
barcode per cell and (b) a single-cell phenotype (transcriptome) per cell.

**Procedure.**
1. Reconstruct the **barcode tree** from the recording array alone. This is the
   ground-truth genealogy; tag `barcode_tree`.
2. Build the **phenotype metric** from transcriptomes alone; tag
   `phenotype_metric`. These two objects share no features.
3. M1: \(\beta\) = distinguishable-lineage growth per division from the barcode
   tree.
4. M2, M3: \(c\) and \(h_{\mathrm{pack}}\) from the phenotype embedding.
5. M4: \(\delta\) of the phenotype metric scored against the barcode tree.
6. M5: \(\eta\) (tags differ, so the meter computes it).
7. Repeat on \(\ge 3\) independent reconstructions (individuals / organoids).

**Decision rule.**
- *Predict:* \(\eta\le 1\) at every depth; if description length is under
  selection, \(\eta\) is high (\(\ge 0.8\)) and stable (across-reconstruction
  SD \(< 0.1\)).
- *Kill (bound):* independently measured \(\beta > c\,h_{\mathrm{pack}}\) at
  fixed \(\varepsilon\) → challenges the bound, so re-run E1 first; if the meter
  is certified, this is a genuine anomaly.
- *Kill (saturation):* \(\eta\ll 1\) systematically → saturation fails for this
  class; the bound stands.

**Why it matters.** Every prior circularity diagnostic in this program traces
to a sequence-inferred tree used on both sides. E3 deletes that shared input.

---

## E4 — Mutation-rate intervention (confound-isolated)

**Question.** Is saturation a *response* to a named knob, not a correlation?
Highest-ranked experiment.

**Inputs.** A serial-passage population with a titratable mutation rate:
mutator alleles (e.g. *mutT/mutS*), a mutagen dose series, or an
error-tunable polymerase in directed evolution. A phenotype or k-mer metric
that does not take mutation rate as input.

**The confound problem (this is the load-bearing addition).** A mutation-rate
knob does not move only \(\beta\). It also shifts the mutational *spectrum*
(transition/transversion bias), drags effective population size \(N_e\), and
changes selection efficiency. A null result is only interpretable if these are
separated from \(\beta\).

**Procedure.**
1. **Assay, don't assume, \(\beta\).** Measure realized substitution rate per
   generation directly (mutation-accumulation lines or deep sequencing), not
   the nominal knob setting. Tag `clock`.
2. **Spectrum control.** Record the transition/transversion and base-composition
   spectrum per condition. Include at least one pair of conditions with
   *matched spectrum but different \(\beta\)* (e.g. two mutagen doses of the
   same agent) and one pair with *matched \(\beta\) but different spectrum*
   (different agents titrated to equal rate).
3. **\(N_e\) control.** Hold census population size and bottleneck schedule
   fixed across conditions; report an \(N_e\) estimate (e.g. from neutral-marker
   variance) and include it as a covariate.
4. Sample at matched generation counts across \(\ge 4\) mutation-rate levels
   spanning \(\ge\)4-fold in \(\beta\).
5. M2, M3, M4 from the representation metric; M5 for \(\eta\).
6. Regress \(c\,h_{\mathrm{pack}}\) on assayed \(\beta\), with spectrum and
   \(N_e\) as covariates.

**Decision rule.**
- *Predict:* \(\partial(c\,h_{\mathrm{pack}})/\partial\beta \approx 1\) in the
  saturated regime, holding spectrum and \(N_e\) fixed; \(\eta\le 1\) always;
  \(\delta\) invariant under the knob (orthogonality control).
- *Kill (saturation-drive):* \(c\,h_{\mathrm{pack}}\) flat under a \(\ge\)4-fold
  \(\beta\) change *with spectrum and \(N_e\) controlled* → the system is not
  driven toward saturation.
- *Confounded, not killed:* if the response appears only when spectrum or
  \(N_e\) also move, the result is inconclusive and the isolation must be
  tightened before any claim.
- *Kill (instrument/bound):* \(\eta>1\) at certified estimators.
- *Kill (orthogonality):* \(\delta\) tracks \(\beta\).

**Near-term proxy.** Time-stamped serially sampled viruses (influenza,
SARS-CoV-2) give \(\beta\) and \(c\) in one physical clock without fossils.
Not an intervention and not spectrum-controlled — run it first as a
**rehearsal of the pipeline**, not as evidence for saturation-drive.

---

## E5 — Description-length pressure in a trained hierarchy (adversarial)

**Question.** Is saturation cross-domain — does a knob that *is* description
length, in a non-biological host, drive \(\eta\to 1\)?

**The false-positive problem (this is the load-bearing addition).** If the same
hand designs a hierarchical generator and a hierarchical loss, high \(\eta\) can
be baked in — you would measure that generator and loss agree, not that the law
holds. E5 is decisive only with an adversarial design and a real negative
control, promoted here to co-equal arms.

**Inputs.** A synthetic generator with tunable alphabet entropy \(h\); a family
of training objectives; a representation whose geometry is measured by M2–M3.

**Arms (all run, pre-registered together).**
1. **Positive arm.** Hierarchical data + a lossless/near-lossless objective
   (rate-distortion, \(\beta\)-VAE, hierarchical contrastive).
2. **Loss-strength sweep.** The same, across rate-distortion parameter values,
   predicting an *ordered* \(\eta\) plateau (monotone in the rate knob).
3. **Negative control — no pressure.** The same data + an objective with **no
   description-length term** (pure reconstruction at unconstrained capacity, or
   a shuffled-target loss). Must **fail** to saturate.
4. **Negative control — no hierarchy.** Non-hierarchical (flat / cyclic) data +
   the lossless objective. Must **fail** to produce exponential retained
   novelty at finite \(c\).
5. **Firewall.** \(\beta\) is read from the *generator* (tag `generator`), never
   from the network; M2–M3 from the representation (tag `representation_metric`).

**Procedure.** Train each arm to convergence across \(\ge 3\) seeds; record
\(\eta(t)\) per epoch; compare plateaus across arms.

**Decision rule.**
- *Predict:* arm 1 \(\eta\to\) near 1; arm 2 plateaus ordered by the rate knob;
  arm 3 and arm 4 plateau well below 1 (the controls fail, as required);
  Euclidean-constrained architectures cannot host arm-1 novelty at finite \(c\)
  (Corollary 4.3: \(\beta\to0\), \(c\to\infty\), or collapse).
- *Kill (saturation-drive):* arm 1 lossless \(\eta\) well below 1 after
  finite-size control.
- *Kill (baked-in artifact):* the negative controls (arms 3–4) also saturate —
  then arm 1's high \(\eta\) is an artifact of the design, not evidence, and the
  experiment proves nothing until the controls are made to fail.
- *Kill (polynomial exclusion):* a Euclidean host successfully carries retained
  exponential novelty at finite \(c\).

---

## E6 — Radius identification (the \(c\)-axis)

**Question.** Is radius accumulated information rather than elapsed time?

**Inputs.** A clade sample where process time, sequence divergence, and
functional complexity are known to disagree: bradytelic ("living fossil")
lineages (deep time, low change) and rapid-innovation lineages (modest time,
high functional gain).

**Procedure.**
1. Process time from fossils or a clock that is **not** the embedding radius.
2. Sequence divergence from an alignment that is **not** the embedding.
3. Functional complexity from an independent measure (gene-content, proteome
   functional-category count, morphological character count).
4. Embedded radius \(r\) via M2 from a representation taking none of the above
   as input.
5. Partial correlations of \(r\) with each of time, divergence, complexity,
   controlling for the other two.

**Decision rule.**
- *Predict:* \(r\) tracks functional complexity (partial \(r>0.5\)) and does
  **not** track elapsed time (partial \(r\approx 0\)) when the three disagree.
- *Kill:* \(r\) tracks elapsed time rather than complexity → the identification
  of \(c\) with an information rate fails; the bound is untouched and the CCS
  radial axis stays advisory.

---

## E7 — Reticulation intervention (\(\delta\perp h_{\mathrm{pack}}\) in vivo)

**Question.** Do volume and thinness stay orthogonal on living cells?

**Inputs.** A microbial population with titratable horizontal transfer:
conjugative plasmids at donor-density series, or a phage-transduction gradient.
A marked mobile element for an independent transfer count.

**Procedure.**
1. Titrate transfer rate across \(\ge 4\) levels; hold the substitution process
   as fixed as the apparatus allows (report it).
2. Transfer rate from marker-acquisition counts (tag `transfer_marker`), not
   from the phenotype metric.
3. M3 (\(h_{\mathrm{pack}}\)) and M4 (\(\delta\)) from the chromosome k-mer /
   phenotype metric.
4. Regress each on transfer rate.

**Decision rule.**
- *Predict:* \(\delta\) increases with transfer rate; \(h_{\mathrm{pack}}\)
  approximately invariant (the synthetic \(\delta\perp h_{\mathrm{pack}}\) split
  certified in E1's M4 run, reproduced in vivo).
- *Kill:* \(\delta\) blind to a \(\ge\)4-fold transfer change, or
  \(h_{\mathrm{pack}}\) tracking transfer as if it were branching → operational
  orthogonality fails for this class.

---

## E8 — Boundary mapping (where the law stops)

**Question.** Do the premises do real work — is there a hierarchy that creates
and retains yet provably falls outside the theory?

**Rationale.** The strongest evidence for a tier-1 law is a clean boundary:
finding *where* it stops applying shows the premises are load-bearing rather
than always-satisfiable. Every other experiment probes systems assumed inside
the theory; E8 hunts the edge.

**Candidate boundary systems.**
1. **Overwriting memory.** A process that changes state fast but does **not**
   retain history (a well-mixed chemostat at equilibrium, a Markov chain with
   short mixing time). Premise 1 (retention) fails → \(\beta\) should be \(\sim
   0\) despite rapid instantaneous change.
2. **Collapsing resolution.** A hierarchy whose distinguishability
   \(\varepsilon(R)\) shrinks with depth (measurement noise growing with time).
   Premise 2 fails → the naive bound should appear violated until
   \(\varepsilon(R)\) is modelled, then restored.
3. **Euclidean-hostable hierarchy.** A genuinely branching but *shallow* or
   *path-like* process whose tree embeds in low-dimensional Euclidean space
   (Corollary 4.3 boundary). \(h_{\mathrm{pack}}\) should be polynomial, and
   \(\beta\) forced to \(0\) at finite \(c\).

**Procedure.** For each candidate, measure the full vector and locate which
premise fails. Show that the theory's *predicted* behavior at that boundary
(zero \(\beta\), or restored bound after modelling \(\varepsilon(R)\), or
polynomial \(h_{\mathrm{pack}}\)) is what occurs.

**Decision rule.**
- *Predict:* at each boundary the named premise fails and the theory's
  degenerate prediction holds (e.g. overwriting → \(\beta\approx0\); resolution
  collapse → apparent violation that the \(\varepsilon(R)\) correction removes).
- *Kill (scope inflation):* a system that satisfies all premises yet violates
  \(\beta\le c\,h_{\mathrm{pack}}\) at certified estimators → challenges the
  bound. A system that fails a premise yet the naive \(\eta\) still behaves as
  if saturated → the premises are not doing the work the theory claims, and the
  applicability story needs revision.

---

## Ranked reading

| Rank | Experiment | What it decides | Class |
|---|---|---|---|
| 1 | E4 mutation-rate intervention | Is saturation a response? | intervention |
| 2 | E5 trained hierarchy (adversarial) | Is saturation cross-domain? | intervention |
| 3 | E3 barcoded lineages | Is \(\eta\) real with a given tree? | ground-truth observation |
| 4 | E7 reticulation intervention | Is \(\delta\perp h_{\mathrm{pack}}\) in vivo? | intervention |
| 5 | E8 boundary mapping | Do the premises do real work? | intervention / observation |
| 6 | E6 radius identification | Is \(c\) information? | pre-registered prediction |
| 7 | E2 numerical achievability | Which half of Conjecture 4.4 to prove? | numerical |
| 8 | E1 meter certification | Are the instruments legal? | calibration |

E1 is last in the table and first in time. Nothing above it is interpretable
until the meter recovers synthetic truth.

## What would constitute firm theoretical ground

Not any one of these. The bound is already firm. Firm *applicability* is:

- E1 passed (legal instruments). **Current state:** M2 and M4 certified; M3
  certified only for the exponential/polynomial dichotomy via the growth-class
  gate, not yet as a magnitude estimator — so \(\eta\) to \(\pm20\%\) is not yet
  interpretable and this precondition is not fully met;
- E3 showing \(\eta\le 1\) with a given tree (premises instantiated, bound
  respected in a living hierarchy);
- E4 showing capacity tracking \(\beta\) with spectrum and \(N_e\) controlled
  (saturation as a response, not a correlation);
- E8 showing a clean boundary (the premises are load-bearing).

Firm *generality* is E5 in addition, with its negative controls failing as
required: the same response under a knob that is description length and is not
biology.

Firm *relational geometry* is a proof or a clean numerical settlement of
Conjecture 4.4 (E2 as reconnaissance).

A single failed intervention does not unwind the theorems. It reclassifies
saturation from a law-like regularity to a domain fact — a legitimate and
fully specified outcome of this protocol. The ceiling is explicit: even a clean
sweep yields "the bound holds, the premises are broadly instantiated, and
selected or trained systems fill the budget in the measured classes." That is
second-law-firm — a statement about a class of systems, never a universal
theorem about nature.
