# The throughline — limit, ladder, classifier, chart

This is the conceptual front door to Active Geometry. It is deliberately
smaller than the full proof spine and the experiment registry.

## One sentence

A process that retains distinguishable histories needs enough geometric room
to address them; a finite-alphabet sequence metric can preserve endpoint
distinguishability while losing genealogical resolution, whereas weighted
tree codes in a hyperbolic host can approach block capacity with no
exponential-order relational tax.

That sentence is the program. It has four load-bearing pieces:

1. **Limit:** retained novelty cannot outgrow addressable room.
2. **Ladder:** distinguishable endpoints and preserved genealogy are different
   capacities.
3. **Classifier:** quartets decide whether a measured metric still carries a
   tree.
4. **Chart:** polar \(\mathbb H^2\) is the minimal-dimensional hyperbolic
   candidate in the stated embedding class.

Saturation asks whether nature fills this host. It is a detachable empirical
question, not a fifth piece.

## The map

```text
                         retained genealogy
                                 |
                 +---------------+---------------+
                 |                               |
        finite-alphabet record             polar H² chart
        (sequence distance)                (encoder: depth + angle)
                 |                               |
       endpoints remain distinct          supremum relational
       while quartets lose                capacity equals block
       resolvability at depth             capacity (weighted clock)
                 |                               |
             BALLOON                     CCS CONSTRUCTION
        addresses without tree        candidate tree-preserving host
```

The left branch is why packing alone is insufficient. The right branch is why
the coordinate-system construction belongs to the theory. They are not two
programs.

## 1. Limit — how much can be remembered

For retained-history growth \(\beta\), radial rate \(c\), and host packing
entropy \(h_{\mathrm{pack}}\),

\[
\boxed{\beta\le c\,h_{\mathrm{pack}}.}
\]

This is curvature-free. It is a counting theorem: an
\(\varepsilon\)-faithful finite-rate representation cannot retain more
histories than its metric balls can separate.

For a proper host \(M\), the ceiling is operational, not merely a converse:

\[
\boxed{
C_{\mathrm{block}}(M,o;c,\varepsilon)
=c\,h_{\mathrm{pack}}(M,o;\varepsilon).
}
\]

Each finite-radius block optimum is attained. Codes may be redesigned at every
depth; the identity does not establish one nested, causal, or relational
capacity-achieving code. It completes the theory of endpoint addressability.

## 2. Ladder — addresses are not genealogy

Block codes preserve distinguishability. Relational codes also preserve the
source hierarchy. Their capacities satisfy

\[
C_{\mathrm{rel}}\le C_{\mathrm{causal}}
\le C_{\mathrm{persistent}}\le C_{\mathrm{block}}.
\]

For a process admissible under the specified relational class, so that
\(\beta\le C_{\mathrm{rel}}\), the slack has an exact accounting:

\[
C_{\mathrm{block}}-\beta
=
\underbrace{(C_{\mathrm{block}}-C_{\mathrm{rel}})}_{\Gamma:
\text{ relational tax}}
+
\underbrace{(C_{\mathrm{rel}}-\beta)}_{\Delta_{\mathrm{use}}:
\text{ utilization slack}}.
\]

A **balloon** has many distinguishable addresses but little relational
capacity. It is the counterexample that prevents “sequence diversity,”
“packing,” and “genealogy” from being treated as synonyms.

In a weighted hyperbolic host the separation closes:

\[
\boxed{
C_{\mathrm{rel}}^{\mathrm{wt}}(\mathbb H_\kappa^n,c;\varepsilon)
=c(n-1)\sqrt\kappa
=C_{\mathrm{block}}(\mathbb H_\kappa^n,o;c,\varepsilon).
}
\]

The supremum over weighted-clock, quasi-isometrically embedded tree codes
equals block capacity. Thus the host's relational tax is zero at exponential
order. This does not establish boundary attainment, equal-edge achievability,
or zero tax for every prescribed genealogy, and it must not be transferred to
an arbitrary sequence metric.

## 3. Classifier — quartets are the right angle

For every quartet \(a,b,c,d\), form the three pair sums

\[
d(a,b)+d(c,d),\quad d(a,c)+d(b,d),\quad d(a,d)+d(b,c).
\]

An exact finite metric is an additive tree metric exactly when the two largest
sums are equal. For noisy or quasi-isometric representations, quartet defect
is a diagnostic requiring a declared tolerance and uncertainty model. It asks
a question orthogonal to packing:

- packing asks **how many addresses fit**;
- quartets ask **whether their relations still form a tree**.

Quartets classify; they do not calibrate \(\kappa\), capacity, or saturation.
They must be evaluated on a representation metric not manufactured from the
same inferred tree. An unresolved quartet must also be distinguished from a
confidently wrong quartet: channel exhaustion can destroy statistical
resolvability without making the pair sums impossible to compute.

## 4. Chart — the coordinate-system construction

Within connected smooth Riemannian hosts, a genuinely branching tree cannot
embed isometrically in dimension one. In the stated embedding class,
\(\mathbb H^2\) is therefore the minimal-dimensional hyperbolic candidate;
path trees remain one-dimensional. In the polar encoder, interpreting radius
as process depth and angle as divergence is a modeling choice, not a
consequence of the capacity theorem.

This is the durable content of the older canonical-coordinate-system
exploration, retained in the
[`minimal_encoder`](../experiments/minimal_encoder/README.md):

- inhabit \(\mathbb H^2\) rather than infer \(n=2\) from the state equation;
- freeze \(\kappa\), because contrastive temperature and curvature are
  non-identifiable;
- train topology from quartets and radius from a separate depth signal;
- ask whether coordinates agree across seeds up to the expected global
  \(O(2)\) gauge (rotation and reflection), unless orientation is fixed.

Procrustes stability supports reproducibility **within the imposed model**. It
does not compare host classes: \(\mathbb H^2\), curvature, taxonomy
supervision, and the radial target are supplied by construction. It is not
evidence that \(\kappa=5/4\), that radius is accumulated information, or that
the biosphere fills an atlas.

## What the balloon simulation adds

A reported ground-truth simulation on Yule trees supplies a finite-sample
analogue of the ladder splitting. Its witness statistics are not estimates of
the formal asymptotic capacities:

- under JC69, the block witness count approaches \(255/256\) while quartet
  resolvability collapses with depth;
- among still-resolvable quartets, topology accuracy remains about
  \(0.78\)–\(0.80\), so declining resolvability dominates the reported depth
  trend (conditional accuracy alone does not exclude error among resolved
  quartets);
- an event-matched infinite-sites control remains at least \(0.998\)
  accuracy, identifying site reuse and multiple hits as the mechanism;
- the observed resolution boundary follows \(D^*\propto\log L\)
  (\(R^2=0.9962\)) within the reported Yule/JC69 sweep.

This result animates the balloon already present in the mathematics. It does
not establish a regime in a real clade, and it remains a reported local result
until its code and outputs are part of the reproducible repository.

It also separates two uses of “saturation” at the level of the reported finite
witnesses:

- **molecular-channel saturation** reduces relational resolvability while
  endpoint distinguishability remains high;
- **host-utilization saturation** would make relational utilization
  \(\eta_{\mathrm{rel}}:=\beta/C_{\mathrm{rel}}\) approach one.

Block efficiency remains
\(\eta_{\mathrm{block}}:=\beta/C_{\mathrm{block}}\). The two efficiencies
coincide only in the weighted hyperbolic class where the capacities are equal.
Molecular and host-utilization saturation are different directions on the
capacity ladder.

## The detachable question

If an attainable host-capacity cost is strictly increasing and lossless
addressability is required, its least-capacity minimizer lies on

\[
h_{\mathrm{pack}}^*=\beta/c.
\]

This variational result makes \(\eta_{\mathrm{block}}=1\) the optimum
**conditional on that objective**. Whether evolution supplies such a cost is
empirical. The state equation and radial concentration belong here,
downstream of both a verified hyperbolic host and measured near-capacity
utilization. Nothing above depends on nature saturating. The sequel contract
is [`IIB_CONTRACT.md`](IIB_CONTRACT.md).

## What is not part of the throughline

- finite-sample growth-class refusal is an instrument guardrail, not “Layer
  0”;
- radial concentration is a consequence of host class plus near-capacity
  utilization, not a foundation;
- \(\kappa=5/4\) from genetic-code entropy is not a theorem;
- \(n=2\) is not back-solved from curvature;
- a coordinate chart is not a filled atlas or an information clock;
- selection does not imply saturation until a biological capacity cost is
  identified.

## Stop here

The four-piece spine is closed at the block and weighted-clock supremum level.
Inside the Heintze class, Theorem 7.1 converts the old genericity conjecture
into A3 \(\Rightarrow\) real hyperbolic space. A3 itself is asserted, not
measured. Equal-edge synchronization, host-class selection (E9), and
biological saturation remain open. The next work is measurement, not another
layer:

1. land and reproduce the balloon simulation;
2. test block witness count and quartet resolvability on an independently
   referenced deep alignment;
3. run the matched-capacity Euclidean-versus-hyperbolic intervention (E9);
4. test the saturation drive separately (E5) only if the Layer IIb question is
   being pursued — and only after the embedder recovers exponential fan-out
   ([`experiments/E5_DESIGN_FLAW.md`](../experiments/E5_DESIGN_FLAW.md)).

The detailed proof and claim status remain in
[`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md) and
[`CLAIMS.md`](CLAIMS.md). This page is only their minimal dependency graph.
