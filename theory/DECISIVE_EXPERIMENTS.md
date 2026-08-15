# Decisive Experiments for the Addressability Limit

## What an experiment can and cannot establish

The theorems are not at stake. The packing converse and the block identity are
counting facts, machine-checked, and no measurement can strengthen or damage
them. Experiments decide the three questions mathematics cannot, **one per
layer** (see [`CLAIMS.md`](CLAIMS.md) for the layer definitions):

1. **Applicability (Layer I).** Do natural systems instantiate the premises —
   retained distinguishable histories, fixed resolution, finite radial rate,
   and *exponential* host growth?
2. **Host class (Layer IIa).** Is the host hyperbolic/tree-like, and is that
   forced — does a hyperbolic host beat a Euclidean one *at matched capacity*?
   Isotropy is currently asserted, never measured; this is where it gets tested.
3. **Saturation (Layer IIb).** Where the host is fixed, does the process fill
   its budget (\(\eta\to 1\)), and is it *driven* there?

Decisiveness is ranked: **interventions** (turn a knob the theory names,
predict the response function) beat **pre-registered predictions** (state the
number before measuring) beat **calibrated observations** (fit after the
fact).

This document is written to be executed. Every experiment names its layer (the
canonical assignment is the ranked table at the end), its inputs, a numbered
procedure, the estimator it calls (from the shared library below), and a
**decision rule** with explicit numeric thresholds — a predicted outcome and a
kill line. A protocol without a kill line is not an experiment.

The set is nine experiments across the three layers. A prior version tested
**Layer IIb saturation almost exclusively** while the better-supported
biological claim is **Layer IIa host class** — a misallocation. E9
(matched-capacity Euclidean vs hyperbolic), pre-registered elsewhere as
*E-alpha* and never run, is the decisive IIa test and was missing; it is added
here and ranked accordingly. The allocation is audited in the table at the end.

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

**Isotropy is a premise, not a measurement.** The Layer IIa hypothesis that the
host is isotropic hyperbolic appears in the meter only as the flag
`--assume-isotropic-hyperbolic`. No experiment below *asserts* it; E9 is the
one that *tests* it. Any curvature magnitude reported without E9 having passed
is conditional on an unverified switch, and must be labelled so.

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

## E2 — Equal-edge endpoint diagnostic

**Layer IIa refinement, not a test of host capacity.** Theorem 4.4 now proves
that weighted/radial relational capacity of \(\mathbb H_\kappa^n\) equals
\(c(n-1)\sqrt\kappa\). E2 asks a stronger synchronization question: what if
every source-tree edge is forced to consume exactly one clock unit, the budget
is exactly \(cR\) at every finite depth, and no startup slack is allowed?

**Question.** Is the exact boundary rate attained in that equal-edge subclass,
or does synchronization force endpoint distortion?

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
   the former unit-edge conditions of Conjecture 4.4 still hold.
3. Plot \(\hat\beta(R)\) against the bound \(c(n-1)\sqrt\kappa\), and
   \((D,K)\) against \(R\).

**Decision rule.**
- *Endpoint attained:* \(\hat\beta(R)\to c(n-1)\sqrt\kappa\) with
  \((D,K)\) bounded in \(R\).
- *Endpoint obstructed:* \(D(R)\to\infty\) at the exact boundary rate.
- Neither outcome changes Theorem 4.4: capacity is a supremum and the weighted
  coding theorem is already proved.

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

This is the **endpoint-obstructed** signature for this layout, not a host
capacity gap. Theorem 4.4 supplies the decisive subcritical construction:
Bishop--Jones free semigroup trees with local generator-dependent edge
durations have critical exponents arbitrarily close to the ambient exponent
and quasi-isometric orbit maps. The equal-edge restriction is exactly what
E2 adds and what the theorem does not need.

A prior fixed-cone implementation reached an artifactual
\(\varepsilon\)-decay verdict; that two implementers of the same prose diverged
is why step 1 pins the construction.

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

**Layer IIa + IIb consequence, not an unconditional prediction.** Theorem 5.1
of [`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md) proves:
*if* the host is hyperbolic (IIa) and the code operates near capacity (IIb),
then almost all retained histories use nearly all of the radial budget
available at **their own weighted clock time**. At exact capacity,

\[
\frac{d(o,f(v))}{c\tau(v)}\to1
\quad\text{in probability under uniform counting on clock balls}
\]

along a capacity-realizing sequence (or eventually when the growth limit
exists). With additive capacity deficit
\(\Delta_{\rm cap}=ch-\beta=ch(1-\eta)>0\), the exact statement is: for every
margin \(m>0\), the fraction with
\[
\frac{d}{c\tau}\le
1-\frac{\Delta_{\rm cap}}{hc}-m
=\eta-m
\]
vanishes along the applicable sequence. Existing independent saturation tests
fail, so the exact-capacity prediction is **not currently licensed for biology
without a separate \(\eta\) measurement**.

**Question.** Conditional on a hyperbolic host and measured near-capacity
operation, does biological radius show the theorem's clock-relative
concentration?

**Inputs.** A hierarchy with an independently calibrated **weighted generative
clock \(\tau\)**, a representation radius, and enough retained histories at
successive clock cutoffs to estimate threshold fractions. A barcoded lineage
system (E3) is preferable. Fossil elapsed time, sequence divergence, or a
complexity score is not automatically the theorem's \(\tau\); substituting one
requires an independent calibration.

**Procedure.**
1. Define \(\tau\) independently of the embedding and verify that clock balls
   \(T_R^\tau\) are finite and sampled without lineage bias (or reweight to the
   uniform counting measure).
2. Estimate embedded radius \(d(o,f(v))\), radial rate \(c\), offset \(A_0\),
   and host exponent \(h\) from independent inputs.
3. Independently verify the theorem's premises: fixed resolution,
   \(\varepsilon\)-separation, quasi-isometry/relational fidelity, the radial
   budget, and hyperbolic host class (E9).
4. Independently estimate growth \(\beta\), then
   \(\Delta_{\rm cap}=ch-\beta\) and efficiency \(\eta=\beta/(ch)\).
   Without both, label the run exploratory.
5. For pre-registered \(\delta\) values and successive \(R\), compute
   \[
   q_\delta(R)=
   \frac{\#\{v\in T_R^\tau:
   d(o,f(v))\le(1-\delta)c\tau(v)+A_0\}}
   {|T_R^\tau|}.
   \]
   Test its decay against \(h\delta c-\Delta_{\rm cap}\).
6. Only as a secondary biological interpretation, compare radius with elapsed
   time, sequence divergence, and functional complexity where they disagree.

**Decision rule.**
- *Necessary-condition prediction:* when E9 supports the hyperbolic host and
  independently measured \(\eta\approx1\), \(q_\delta(R)\) decays
  exponentially; at exact capacity \(r/(c\tau)\) concentrates near 1.
- *Positive result:* supports a necessary consequence of near-capacity
  hyperbolic coding, but does **not** prove saturation.
- *Negative result:* challenges the biological application only after
  separation, quasi-isometry, radial budget, unbiased counting, asymptotic
  regime, host class, and near-capacity premises are independently verified.
  Otherwise it is a premise failure, not a theorem failure. Use E9 and the
  independent \(\eta\) measurement to separate causes. The bound is untouched
  and the CCS radial axis stays advisory.

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

## E9 — Matched-capacity Euclidean vs hyperbolic (the realization test)

**Layer IIa — host class. This is the decisive realization test, and the one
the program most needs.** Pre-registered elsewhere as *E-alpha* ("Matched-
Capacity Euclidean vs Hyperbolic," deferred to a GPU phase) and **never run**;
absent from the prior protocol. It is distinct from E5's Euclidean control,
which tests Layer I polynomial exclusion (*can* a Euclidean host carry
exponential novelty at all), not *whether hyperbolic is forced*.

**Question.** Given the same data and a matched capacity budget, does a
hyperbolic host represent tree-structured biology with materially higher
relational fidelity than a Euclidean host? If yes, the hyperbolic host class is
forced by the data's structure (Conjecture 7.1's empirical shadow). If they
tie, the hyperbolic realization is *not* forced and Layer IIa is unsupported —
however well saturation (IIb) fares.

**Inputs.** A tree-structured dataset with a ground-truth or barcode genealogy
(share E3's data where possible); two embedding hosts, \(\mathbb E^m\) and
\(\mathbb H^n\); a capacity-matching protocol.

**Capacity matching (the load-bearing design).** "Matched capacity" must be
matched *packing capacity*, not matched parameter count. Fix resolution
\(\varepsilon\) and radius budget, then size the two hosts so their certified
block capacities \(c\,h_{\mathrm{pack}}\) (M3 growth-class-gated) agree within
error. Matching parameters instead of capacity would confound the test with
representational budget.

**Procedure.**
1. Match capacities as above; verify with the E1-certified M3 gate that both
   hosts are in the same growth class at the matched budget.
2. Embed the same genealogy into each host under identical relational
   objectives and seeds.
3. Measure relational fidelity in each: distortion \((D,K)\), tree-defect M4
   against the ground-truth tree, and held-out lineage-placement accuracy.
4. Repeat across \(\ge 3\) seeds and \(\ge 2\) dimensions per host.

**Decision rule.**
- *Predict (IIa holds):* hyperbolic achieves lower distortion and defect at
  matched capacity, with a gap that widens with tree depth; Euclidean fidelity
  degrades as depth grows (the crowding the packing bound predicts).
- *Kill (IIa fails):* Euclidean matches hyperbolic fidelity at matched capacity
  across depths → the hyperbolic host class is not forced; the realization
  claim, and with it the biological use of curvature, is unsupported. Neither
  the Layer I bound nor the block identity is touched.

**Status.** Not run; needs GPU-scale embedding. It is the highest-value
unrun experiment in the program, because the biological claim is IIa and E9 is
the only direct IIa intervention.

---

## The empirical bridge: the growth-class × tree-defect phase diagram

The two certified instruments (the growth-class gate on M3, the defect meter
M4) yield a two-axis classifier that needs **no curvature magnitude** and no
uncertified \(\eta\). It is the honest empirical figure for a real system, and
it maps a dataset onto the two-layer theory directly.

| | tree-like: \(\delta\approx 0\) | reticulate: \(\delta>0\) |
|---|---|---|
| **exponential capacity** | hyperbolic-tree host candidate (Layer II applies) | exponential network / mixed host |
| **polynomial capacity** | path-like or non-expanding hierarchy | networked but capacity-limited |

Reading it asks exactly the two questions the layers separate:

1. *Growth class* (M3 gate): does the metric supply exponential room? This is
   Layer I — Corollary 4.3, mechanically enforceable today.
2. *Tree defect* (M4): is the relational structure tree-like, and if not, how
   far (the measured reticulation residual)?

Only the upper-left cell motivates a hyperbolic realization, and even there it
does **not** fix a curvature magnitude — that needs a certified M3 magnitude
estimator and the independence firewall. The disciplined progression for any
substrate is therefore: **growth class → tree defect → constrained capacity →
utilization (saturation) → optional curvature realization.** Leading with
\(\kappa\) inverts this order and is what produced the retired v2 claims.

---

## Ranked reading

Ranked by decisiveness *and* by layer, so the realization layer (IIa) — where
the biological claim actually lives — is no longer buried.

| Rank | Experiment | Layer | What it decides | Class |
|---|---|---|---|---|
| 1 | E9 matched-capacity Euclidean vs hyperbolic | IIa | Is the hyperbolic host *forced*? | intervention (unrun) |
| 2 | E4 mutation-rate intervention | IIb | Is saturation a response? | intervention |
| 3 | E3 barcoded lineages | IIb | Is \(\eta\) real with a given tree? | ground-truth observation |
| 4 | E7 reticulation intervention | IIa | Is \(\delta\perp h_{\mathrm{pack}}\) in vivo? | intervention |
| 5 | E2 equal-edge endpoint | IIa refinement | Does exact synchronization obstruct the endpoint? | numerical |
| 6 | E5 trained hierarchy (adversarial) | IIb | Is saturation cross-domain? | intervention |
| 7 | E8 boundary mapping | I | Do the premises do real work? | intervention / observation |
| 8 | E6 radius concentration | IIa+IIb | Does clock-relative radius concentrate as the converse predicts? | conditional prediction |
| 9 | E1 meter certification | I | Are the instruments legal? | calibration |

E1 is last in the table and first in time. Nothing above it is interpretable
until the meter recovers synthetic truth.

### Allocation audit (the misallocation, fixed)

The claim best supported by existing biology is the **host class (IIa)** —
seed-stable embeddings, curvature as a fixed design parameter, tree-defect near
zero, \(n=2\) as embeddability. Yet the prior protocol pointed almost entirely
at **saturation (IIb)**: of E1–E8 only E2 and E7 touched IIa, and E2 sat 7th.
The decisive IIa intervention, E9, was missing. This is now corrected: E9 is
added and ranked first, E2 and E7 are read as IIa, and the layer column makes
the balance visible. The standing imbalance in *evidence* remains and is stated
plainly — IIa is comparatively well-supported and under-tested by intervention;
IIb is heavily tested and fails its independent kill lines (domain
\(r=0.35\), protein per-family \(-0.11\), tree-independent viral \(0.06\)–\(0.19\)).

## What would constitute firm theoretical ground

Not any one of these. The bound is already firm. Firm ground is layer by layer.

**Layer I — applicability.**
- E1 passed (legal instruments). **Current state:** M2 and M4 certified; M3
  certified only for the exponential/polynomial dichotomy via the growth-class
  gate, not yet as a magnitude estimator — so \(\eta\) to \(\pm20\%\) is not yet
  interpretable and this precondition is not fully met;
- E8 showing a clean boundary (the premises are load-bearing).

**Layer IIa — host class (the biological claim).**
- E9 showing hyperbolic beats Euclidean at matched capacity — the isotropy
  premise earned rather than asserted. **This is the single most important
  unrun experiment.** Mathematical support: Theorem 4.4 proves full weighted
  relational capacity in the hyperbolic host. Empirical support: E7
  orthogonality in vivo, the tree-defect and seed-stability evidence. E2 is
  only an equal-edge endpoint refinement.

**Layer IIb — saturation.**
- E3 showing \(\eta\le 1\) with a given tree, and E4 showing capacity tracking
  \(\beta\) with spectrum and \(N_e\) controlled (saturation as a response).
  Current independent tests fail here; relocating saturation to IIb does not
  rescue it.

Firm *generality* is E5 in addition, with its negative controls failing as
required: the same response under a knob that is description length and is not
biology.

Firm *relational capacity* is now mathematical: Theorem 4.4 proves that
weighted/radial genealogy-preserving codes attain the full hyperbolic exponent
as a supremum. E2 concerns only the stronger equal-edge synchronization class.

A single failed intervention does not unwind the theorems. It reclassifies
saturation from a law-like regularity to a domain fact — a legitimate and
fully specified outcome of this protocol. The ceiling is explicit: even a clean
sweep yields "the bound holds, the premises are broadly instantiated, and
selected or trained systems fill the budget in the measured classes." That is
second-law-firm — a statement about a class of systems, never a universal
theorem about nature.
