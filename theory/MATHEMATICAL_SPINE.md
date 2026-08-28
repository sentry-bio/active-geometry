# The Addressability Limit for Information-Generating Hierarchies

## Status and purpose

The compact dependency structure is stated in
[`ADDRESSABILITY_KERNEL.md`](ADDRESSABILITY_KERNEL.md). This document supplies
the full definitions, proof, scope, and measurement consequences.

The theory is organized in **two layers** that must not be conflated.

**Layer I — the universal capacity theory (curvature-free), §§1–5.** Three
independent quantities; a **packing theorem** that forces the capacity
inequality \(\beta\le c\,h_{\mathrm{pack}}\); the exact **block identity**
\(C_{\mathrm{block}}=c\,h_{\mathrm{pack}}\); a **constrained-capacity ladder**
with block as its universal ceiling; and the decomposition of slack into a
relational tax and a utilization term. No curvature, tree, or biology appears.

**Layer II — the curvature realization (where real systems live), §§7–10.** This
layer splits into two sublayers with very different evidential standing.

- **Layer IIa — host class (realization).** *Which* geometry hosts the data:
  the space-form argument and **isotropic realization** (§7), **Theorem 7.1**
  (Heintze isotropy conditional on axiom A3: full \(O(d)\) directional
  symmetry), and \(n=2\) as an embeddability floor. This is the
  **better-supported** biological claim. A3 itself is asserted, not measured.
- **Layer IIb — saturation.** *Whether* a process fills its budget: the
  **capacity-saturation hypothesis** (\(\eta=1\), §5) which, combined with the
  IIa realization, yields the **state equation** (§8). This is the **harder,
  less-supported** claim; every independent test to date sits below its kill
  line.

Only Layer I's inequality is unconditional. The host class is IIa; the state
equation is IIb, and follows only after both saturation and isotropy are
supplied. The four-point tree theorem (§6) is an independent classifier that
bridges to IIa: it can motivate a hyperbolic host class or minimal embedding
dimension, but it does not calibrate capacity or curvature. Sections 11–13
(measurement, falsification, minimal statement) span all layers.

A consequence the split makes sharp: **isotropy is asserted, never measured.**
Corollary 4.3 already excludes polynomial hosts at finite radial rate and
fixed dimension; that is Layer I, not a realization test. E9 of
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md) is a finite-sample
shadow of that corollary at matched packing capacity, distinct from E5's
Euclidean control and not the argument that forces \(\mathbb H^2\). A3
remains the named symmetry axiom for Theorem 7.1.

All logarithms in this document are natural unless explicitly written
`log₂`.

---

# Layer I — The universal capacity theory

## 1. Information-generating histories

Let \(\mathcal H_R\) be the histories retained through generative depth \(R\).
Assume their exponential growth rate exists:

\[
N(R):=|\mathcal H_R|
      =\exp\!\bigl(\beta R+o(R)\bigr),
\qquad
\beta:=
\lim_{R\to\infty}\frac{\log N(R)}{R}.
\]

The quantity \(\beta\) is the production rate of distinguishable retained
history, in nats per generative step. Define the effective retained-information
rate by

\[
h_{\mathrm{eff}}:=\frac{\beta}{\ln2}
\quad\text{bits per generative step},
\]

so that \(\beta=h_{\mathrm{eff}}\ln2\) by definition. A Shannon entropy-rate
estimator is a valid empirical estimator of \(h_{\mathrm{eff}}\) only when a
typical-set or coding argument shows that it measures the exponential growth
of retained distinguishable histories. Alphabet entropy alone does not
establish this.

The word **retained** is load-bearing. A process that overwrites its past need
not have \(N(R)\) growing with \(R\), even if its instantaneous state changes
rapidly.

### Resolution

Fix a distinguishability scale \(\varepsilon>0\). Distinct histories need only
remain distinguishable at this fixed operational resolution. If resolution
instead sharpens with depth, \(\varepsilon=\varepsilon(R)\), its exponential
rate contributes to the bound and must be modeled separately. In a regular
local dimension \(d\), \(\varepsilon(R)\asymp e^{-aR}\) typically adds
\(daR\) to logarithmic packing counts; in a general metric space the correction
is controlled by its small-scale packing dimension.

---

## 2. Representations and their radial rate

Let \((M,d,o)\) be a pointed metric space. A representation

\[
f:\bigcup_R\mathcal H_R\longrightarrow M
\]

is **\(\varepsilon\)-faithful at finite radial rate** when:

1. \(d(f(x),f(y))\ge\varepsilon\) for distinct retained histories \(x,y\);
2. \(f(\mathcal H_R)\subseteq B(o,r(R))\);
3. the radial conversion rate is finite:

   \[
   c:=\limsup_{R\to\infty}\frac{r(R)}R<\infty.
   \]

The units of \(c\) are geometric length per generative step. It converts the
process clock into the host's radial ruler. It cannot be silently omitted when
curvatures measured in different radial units are compared.

The common convention \(c=1\) is the **process-time gauge**: one generative step
is declared to be one radial unit. This is a useful coordinate choice, not a
theorem.

---

## 3. Metric addressing capacity

Let \(P(B(o,\rho),\varepsilon)\) be the supremum of the cardinalities of
\(\varepsilon\)-separated subsets of \(B(o,\rho)\), possibly \(+\infty\).
Define the packing entropy

\[
h_{\mathrm{pack}}(M,o;\varepsilon)
:=
\limsup_{\rho\to\infty}
\frac{\log P(B(o,\rho),\varepsilon)}{\rho}.
\]

This definition needs only a metric. It therefore applies to trees, graphs,
statistical manifolds, learned representations, and Riemannian manifolds.

For a complete Riemannian manifold with uniform positive lower and finite
upper volume bounds for balls of the relevant fixed radii, fixed-resolution
packing entropy agrees with volume entropy at exponential order:

\[
h_{\mathrm{pack}}
=
h_{\mathrm{vol}}
:=
\limsup_{\rho\to\infty}
\frac{\log\operatorname{Vol}B(o,\rho)}{\rho}.
\]

The stated ball-volume bounds make packing counts and ball volume differ only
by multiplicative constants at fixed resolution. We refer to this explicit
condition below rather than using an unspecified notion of bounded geometry.

---

## 4. The addressability theorem

### Theorem 4.1 — Addressability bound

If a hierarchy of growth rate \(\beta\) admits an \(\varepsilon\)-faithful
representation at finite radial rate \(0<c<\infty\), with
\(r(R)\to\infty\), in a host of finite packing entropy
\(h_{\mathrm{pack}}\), then

\[
\boxed{\beta\le c\,h_{\mathrm{pack}}.}
\]

For a Riemannian host satisfying the fixed-ball volume bounds of Section 3,

\[
\boxed{\beta\le c\,h_{\mathrm{vol}}.}
\]

If \(h_{\mathrm{pack}}=+\infty\), the packing inequality is immediate; the
finite-entropy hypothesis is used only for the non-trivial rate comparison.

#### Proof

The represented histories form an \(\varepsilon\)-separated subset of
\(B(o,r(R))\), so

\[
N(R)\le P(B(o,r(R)),\varepsilon).
\]

For every \(\delta>0\), eventually
\(r(R)\le(c+\delta)R\). Monotonicity of packing numbers gives

\[
P(B(o,r(R)),\varepsilon)
\le P(B(o,(c+\delta)R),\varepsilon).
\]

Because \(h_{\mathrm{pack}}\) is a finite limsup, at all sufficiently large
radii

\[
\log P(B(o,\rho),\varepsilon)
\le (h_{\mathrm{pack}}+\delta)\rho
\]

. Hence

\[
\frac{\log N(R)}R
\le
(c+\delta)(h_{\mathrm{pack}}+\delta)+o(1).
\]

Take the upper limit and then let \(\delta\downarrow0\). \(\square\)

The finite-depth packing count and the convergent finite-rate case are
machine-checked as `addressability_limit` in
[`lean/ActiveGeometry/Packing.lean`](lean/ActiveGeometry/Packing.lean), using
Mathlib's exact `Metric.packingNumber`. The limsup statement above is retained
as the more general paper theorem.

### Theorem 4.2 — General block-addressability identity

The converse closes to an identity in full generality when the operational
task is matched exactly to metric packing. Define

\[
A_{\rm block}(\rho,\varepsilon)
:=
\max\left\{
|C|:\ C\subseteq B(o,\rho),\
d(x,y)\ge\varepsilon\ \text{for }x\ne y
\right\}.
\]

For every proper metric host, closed balls are compact and their
fixed-resolution packing numbers are finite and attained. Therefore

\[
\boxed{
A_{\rm block}(\rho,\varepsilon)
=P(B(o,\rho),\varepsilon).
}
\]

The left side is operational: the largest noiseless message set that can be
assigned distinguishable addresses within budget \(\rho\). The right side is
geometric. They are one quantity.

At a linear radial budget \(\rho=cR\), \(c>0\), define

\[
C_{\rm block}(c,\varepsilon)
:=
\limsup_{R\to\infty}
\frac{\log A_{\rm block}(cR,\varepsilon)}R.
\]

Because ball-packing counts are monotone in radius, sampling them at the
linearly spaced radii \(cR\) preserves their radial limsup. Hence

\[
\boxed{
C_{\rm block}(c,\varepsilon)
=c\,h_{\rm pack}(M,o;\varepsilon).
}
\]

This is a genuine operational-geometric capacity identity for arbitrary
proper metric spaces. It requires no tree, hyperbolicity, isotropy,
stationarity, or curvature. The finite-radius existence half is
machine-checked as `exists_optimal_blockCode`; combined with
`card_le_packingCount`, it establishes exact finite-block achievability.

#### Why this does not make every hierarchy capacity-achieving

A block code may be redesigned at every depth and must preserve only terminal
distinguishability. Stronger meanings of “represent a hierarchy” form an
admissibility ladder:

1. **Block:** terminal messages are \(\varepsilon\)-separated.
2. **Persistent:** one nested address family retains prior messages.
3. **Causal:** parent-child addresses also obey a per-step motion or locality
   bound.
4. **Relational:** the representation preserves the hierarchy's pairwise
   metric or topology up to a stated distortion, for example by a uniform
   quasi-isometric embedding.

The corresponding constrained capacities cannot exceed \(C_{\rm block}\).
Packing alone proves achievability only at level 1. Extension theorems,
online-coding arguments, or geometric embedding theorems are needed at the
stronger levels.

Positive packing entropy is insufficient for relational achievability. For a
concrete counterexample, take a ray of hubs and attach to the hub at radius
\(R\) a clique of \(\lceil e^{hR}\rceil\) vertices. With the graph metric this
is a proper space with packing entropy \(h\): exponentially many distinguishable
addresses occur in each bounded-diameter “balloon.” But it cannot contain a
uniform quasi-isometric image of a regular branching tree. Up to radius \(R\)
there are only \(O(R)\) bounded-diameter shells, while a tree level has
exponentially many mutually tree-distant vertices; a quasi-isometric lower
bound permits only a uniformly bounded number of them in each shell.

This example locates the missing hypothesis precisely. Ordinary packing
capacity counts addresses; it does not certify that their mutual geometry
realizes the source's genealogy.

Sarkar-type hyperbolic tree embeddings address the stronger relational task,
but a finite-tree low-distortion result is not by itself an asymptotic
capacity theorem. Closing that theorem requires a uniform family at fixed host
geometry, fixed resolution, controlled radial rate, and distortion that does
not degrade with depth.

Finally, quartet defect \(\delta\) is not the capacity slack
\(\Delta=c\,h_{\rm pack}-\beta\). It classifies deviation from tree geometry.
A reticulate source can use every available address, and a perfect tree can
use very few. Any theorem turning \(\delta\) into a quantitative capacity gap
must add a particular source class, host class, and distortion criterion.

### Theorem 4.4 — Relational capacity of hyperbolic hosts

The problem formerly stated here as Conjecture 4.4 is now settled. The exact
unit-edge formulation was false as quantified: if \(c<\varepsilon\), a child
of the root would have to lie both within distance \(c\) and at least
\(\varepsilon\) from it. The only admissible tree is then the root, so its
capacity is zero rather than \(c(n-1)\sqrt\kappa\).

The error was not geometric; it was clocking. It forced every source-tree edge
to cost exactly one unit while allowing no startup slack, even though the
general addressability theorem permits an independently calibrated process
clock. Give each tree edge \(e\) a positive duration \(a(e)\), let
\(\tau(v)\) be elapsed path time, and use the weighted path metric \(d_\tau\).
A weighted relational code is \(\varepsilon\)-separated, quasi-isometric for
\(d_\tau\), and obeys

\[
d(o,f(v))\le c\,\tau(v)+A_0
\]

for a depth-independent \(A_0\). Define \(C_{\rm rel}^{\rm wt}\) from the
growth of \(\{v:\tau(v)\le R\}\). Then, for every \(\kappa>0\), \(n\ge2\),
\(c>0\), and \(\varepsilon>0\),

\[
\boxed{
C_{\rm rel}^{\rm wt}(\mathbb H_\kappa^n,c;\varepsilon)
=c(n-1)\sqrt\kappa .
}
\]

The upper bound is Theorem 4.1. For the lower bound, choose an
\(\varepsilon\)-separated finite-index subgroup of a cocompact hyperbolic
lattice. Skenderi's 2026 Bishop--Jones theorem supplies free subsemigroups with
critical exponent arbitrarily close to the lattice exponent
\((n-1)\sqrt\kappa\), and its orbit tree is quasi-isometrically embedded.
Assign generator \(s\) the edge duration
\(a_s=d(o,so)/c\). The triangle inequality gives the radial budget exactly,
and Skenderi's partition-sum inequality gives weighted growth arbitrarily
close to \(c(n-1)\sqrt\kappa\).

Thus preserving genealogy has **zero relational tax at exponential order** in
real hyperbolic space. The proof does not require or establish endpoint
attainment. The stronger
equal-edge, exact finite-depth problem remains false for \(c<\varepsilon\) and
unresolved for \(c\ge\varepsilon\); E2 probes that synchronization constraint,
not host capacity.

The complete definitions, counterexample, proof, citation dependencies, and
formal boundary are in
[`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md). The lower
bound uses Aleksander Skenderi, "Free semigroups of large critical exponent,"
*Journal of Topology* (2026), DOI
[10.1112/topo.70087](https://doi.org/10.1112/topo.70087), especially Theorem
3.1, Proposition 3.2(4), Theorem 1.1 (5.1), and equation (6.5).

### Corollary 4.3 — Polynomial-growth exclusion (the trichotomy)

If \(\beta>0\) and \(c<\infty\), then \(h_{\mathrm{pack}}>0\). Consequently, no
polynomial-growth host can represent retained exponential novelty faithfully
at finite radial rate.

The same packing count, with the quantifiers flipped, is the content of the
limit, not a second theorem:

1. **Exponential room.** Retain novelty at finite address cost \(\Rightarrow\)
   \(h_{\mathrm{pack}}>0\).
2. **Exponential addresses.** In a host whose packing grows as a Euclidean
   (or Bishop–Gromov, \(\operatorname{Ric}\ge 0\)) polynomial of degree
   \(n\), faithfulness forces
   \(r(R)\ge\varepsilon\,e^{\beta R/n}(1-o(1))\). Raising \(n\) raises the
   polynomial budget; the exclusion is at *fixed* dimension.
3. **Forgetting.** If the process overwrites rather than accretes,
   \(\beta=0\) and a polynomial host is legal. That exemption is the
   bound's edge, not a counterexample.

Every finite-depth census is a packing number on a compact ball of a proper
host (Theorem 4.2). The displayed rates are the \(R\to\infty\) limit of
those counts; the remainder is \(O(1/R)\). Manning's identity is a
dynamical reading of volume entropy on compact manifolds and is not used.
In exponential-growth hosts the leading packing *rate* is robust in
\(\varepsilon\); resolution enters the intercept.

This is the coordinate-free content of the theory:

> Remembering while creating requires an exponential-growth host — or
> exponentially growing addresses, or forgetting.

Positive packing or volume entropy does **not**, by itself, imply constant
negative sectional curvature or Gromov hyperbolicity. Hyperbolic space is a
canonical realization only after additional symmetry assumptions are imposed.
A matched-parameter Euclidean bake-off (extra \(d\), extra polynomial
budget) is not a test of this corollary. A matched-*packing* comparison at
fixed dimension is its finite-sample shadow (E9), optional, and not
load-bearing for Layer I.

---

## 5. Slack, efficiency, and process saturation

Define the non-negative addressability slack

\[
\Delta:=c\,h_{\mathrm{pack}}-\beta\ge0
\]

and, when the denominator is positive, the efficiency

\[
\boxed{
\eta:=
\frac{\beta}{c\,h_{\mathrm{pack}}}
\le1.
}
\]

The efficiency \(\eta\) is invariant under a change of radial length units
because \(c\) and \(h_{\mathrm{pack}}\) transform inversely. The dimensional
slack \(\Delta\) is expressed in nats per generative step and changes if the
generative-step unit itself is redefined.

- \(\eta=1\): capacity is saturated at exponential order;
- \(\eta<1\): the host has excess asymptotic capacity, or the representation
  pays for distortion, anisotropy, finite-depth effects, or inefficient
  addressing.

### The slack decomposition

The block slack conflates two physically different things once a process is
represented under a relational class. With block capacity
\(B:=c\,h_{\mathrm{pack}}\) and relational capacity \(C_{\rm rel}\le B\),

\[
\boxed{
B-\beta
=
\underbrace{(B-C_{\rm rel})}_{\text{relational tax }\Gamma\ \ge 0}
+
\underbrace{(C_{\rm rel}-\beta)}_{\text{utilization slack }\Delta_{\mathrm{use}}\ \ge 0},}
\]

and the block efficiency factors as

\[
\eta_{\mathrm{block}}
=\frac{\beta}{B}
=\underbrace{\frac{C_{\rm rel}}{B}}_{\text{availability}}\cdot
\underbrace{\frac{\beta}{C_{\rm rel}}}_{\text{utilization}} .
\]

The relational tax \(\Gamma\) is the possible price of preserving genealogy in
a general host or under a stronger synchronization constraint. Theorem 4.4
proves that it is **zero at exponential order in real hyperbolic space** when
the process has a local weighted clock:

\[
C_{\rm rel}^{\rm wt}
=c(n-1)\sqrt\kappa
=B.
\]

Thus in the intended hyperbolic realization, relation preservation itself does
not reduce asymptotic capacity; only endpoint attainment and equal-edge
synchronization may cost. The utilization term remains what "saturation"
means for a particular physical process. In hosts outside Theorem 4.4's class,
an observed \(\eta_{\rm block}<1\) can still mix relational tax and utilization,
so the decomposition remains the correct audit.

### Proposition 5.1 — Exponential saturation

Suppose \(0<c<\infty\), \(r(R)\to\infty\), and the finite limits defining
\(\beta\), \(c\), and the packing entropy along \(r(R)\) exist. If

\[
\log P(B(o,r(R)),\varepsilon)-\log N(R)=o(R),
\]

then

\[
\boxed{\beta=c\,h_{\mathrm{pack}}.}
\]

Conversely, equality of these exponential rates says that the represented
histories fill the available address space at exponential order, though
subexponential inefficiency may remain.

### Theorem 5.2 — Radial concentration converse

Let a weighted relational code in \(\mathbb H_\kappa^n\) have
\(h=(n-1)\sqrt\kappa\), radial rate \(c\), offset \(A_0\), and growth at least
\(ch-\Delta_{\rm cap}\), in the finite form
\(|T_R^\tau|\ge B^{-1}e^{(ch-\Delta_{\rm cap})R}\). Here
\(\Delta_{\rm cap}\) is an additive rate deficit, not the dimensionless
efficiency \(\eta=\beta/(ch)\). For \(0<\delta<1\), let \(E_\delta(R)\) be
histories through time \(R\) satisfying

\[
d(o,f(v))\le(1-\delta)c\tau(v)+A_0.
\]

The smaller-ball packing count gives

\[
\boxed{
\frac{|E_\delta(R)|}{|T_R^\tau|}
\le C e^{-(h\delta c-\Delta_{\rm cap})R}
}
\]

whenever the stated growth lower bound holds, with
\(C=C(\varepsilon,n,\kappa,A_0,B)\). For a limsup growth rate, the estimate
holds along a realizing sequence with an arbitrary exponent margin
\(\xi>0\); if the growth limit exists, it holds eventually. No
shell-regularity assumption is needed.

At exact capacity (\(\Delta_{\rm cap}=0\)),

\[
\frac{d(o,f(v))}{c\tau(v)}\to1
\]

in probability under uniform counting on \(T_R^\tau\), along the realizing
sequence (or eventually when the growth limit exists). For fixed
\(\Delta_{\rm cap}>0\), the quantified conclusion is only that every band
strictly below \(1-\Delta_{\rm cap}/(hc)=\eta\) has vanishing mass. Thus
“radius is accumulated information” is a necessary consequence of
**near-capacity hyperbolic coding**, not an unconditional prediction about
nature. The finite-radius sub-ball count is Lean-checked; the proof,
sharpness-construction sketch, and limsup qualifications are in
[`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md) §5.

### Variational interpretation

Fix \(c\). Consider any admissible-representation objective whose geometric
cost is strictly increasing in an attainable host-capacity variable while
lossless representation imposes \(c\,h_{\mathrm{pack}}\ge\beta\). If a minimizer
exists and the boundary capacity is attainable, its least-capacity solution
lies on the boundary:

\[
\min h_{\mathrm{pack}}
\quad\text{subject to}\quad
c\,h_{\mathrm{pack}}\ge\beta
\quad\Longrightarrow\quad
h_{\mathrm{pack}}^*=\frac{\beta}{c}.
\]

This conclusion does not depend on one specially chosen description-length
formula. What remains empirical is whether a physical system is driven toward
such an economical representation.

The saturation condition is analogous to a reversible Carnot cycle or a
capacity-achieving Shannon code. The inequality is the limit; equality is its
ideal realization.

---

# Bridge — the four-point classifier

## 6. The quartet: classification, not calibration

*This section is neither Layer I nor Layer II. It is an independent
classification theorem that bridges them: it decides whether a metric is
tree-like (and so motivates the hyperbolic host class of Layer II) without
supplying \(\beta\), \(c\), saturation, or curvature magnitude. Its defect
\(\delta\) is the measured residual by which a real reticulate hierarchy
departs from the pure tree class.*

For four points \(a,b,c,d\), define

\[
\begin{aligned}
S_1&=d(a,b)+d(c,d),\\
S_2&=d(a,c)+d(b,d),\\
S_3&=d(a,d)+d(b,c).
\end{aligned}
\]

### Theorem 6.1 — Four-point tree condition

A finite metric is an additive weighted-tree metric exactly when, for every
quartet, the two largest values among \(S_1,S_2,S_3\) are equal.

In the geodesic setting, a metric space is \(0\)-hyperbolic exactly when it is
an \(\mathbb R\)-tree.

These statements are the precise bridge between phylogenetic tree metrics and
Gromov hyperbolicity. They determine whether the qualitative object is a tree.
They do not determine a curvature magnitude.

A genuinely branching tree cannot be embedded isometrically in a
one-dimensional connected Riemannian manifold, while finite trees admit
arbitrarily low-distortion embeddings in the hyperbolic plane after the
curvature scale is chosen appropriately. Thus \(n=2\) is the minimal smooth
hyperbolic ambient dimension for genuinely branching tree metrics within this
embedding class. Path trees remain one-dimensional. The four-point condition
alone therefore does not imply \(n=2\).
It is an embeddability statement, not the minimizer of the curvature
functional.

Reticulation, reassortment, horizontal transfer, or source mixing can violate
the exact quartet condition. Their effects should be reported through
non-zero \(\delta\), increased effective boundary dimension, anisotropy, or a
combination of these—not folded into the tree equality by definition.

---

# Layer II — The curvature realization

*Layer IIa (host class) is §7 (realization, Theorem 7.1) with the
\(n=2\) floor of §6. Layer IIb (saturation) is the equality of §8, conditional
on the coordinate-free saturation of §5. §§9–10 refine scope and dynamics.*

## 7. Isotropic hyperbolic realization

Now add a separate canonicality assumption: the host is complete,
homogeneous, isotropic, simply connected, and unbounded. The space-form
classification and positive entropy select

\[
M=\mathbb H_\kappa^n
\]

with sectional curvature \(-\kappa\), \(\kappa>0\).

Its volume entropy is

\[
h_{\mathrm{vol}}(\mathbb H_\kappa^n)
=(n-1)\sqrt{\kappa}.
\]

The addressability theorem therefore gives the **curvature floor**

\[
\beta\le c(n-1)\sqrt{\kappa},
\]

equivalently,

\[
\boxed{
\kappa\ge
\left(\frac{\beta}{c(n-1)}\right)^2.
}
\]

When capacity is saturated,

\[
\boxed{
\kappa^*=
\left(\frac{\beta}{c(n-1)}\right)^2.
}
\]

With \(\beta=h_{\mathrm{eff}}\ln2\),

\[
\boxed{
\kappa^*=
\left(
\frac{h_{\mathrm{eff}}\ln2}{c(n-1)}
\right)^2.
}
\]

For a genuinely branching tree modeled in its minimal smooth hyperbolic
ambient dimension \(n=2\),

\[
\kappa^*=
\left(\frac{h_{\mathrm{eff}}\ln2}{c}\right)^2.
\]

### Theorem 7.1 — Heintze isotropy (conditional on A3)

Layer I is curvature-free, so it is fair to ask why Layer II is about
hyperbolic space rather than an arbitrary exponential-growth host. Most of
that question is already settled. Geometry alone does not hand over
\(\mathbb H_\kappa^n\) among homogeneous negatively curved hosts. One named
axiom about the generator does.

**Already forced, given isotropy.** Among complete, homogeneous, isotropic,
simply connected, unbounded Riemannian hosts, the space-form classification
leaves constant-curvature models. The sphere is excluded by bounded diameter.
Euclidean space is excluded by polynomial packing (Corollary 4.3). Positive
volume entropy and relational realizability of branching trees therefore
select \(\mathbb H_\kappa^n\). Theorem 4.4 then says that, in this host,
genealogy has zero exponential-order tax as a weighted-clock supremum.

**Geometry will not do the rest.** Drop isotropy, keep homogeneity and
negative curvature. Heintze's classification identifies such manifolds with
solvable groups \(N\rtimes\mathbb R\); the generic member has non-abelian
\(N\) or a non-scalar contraction, hence pinched, anisotropic curvature.
Real hyperbolic space is the symmetric special case (\(N\) abelian,
contraction a scalar). The former wording of Conjecture 7.1 — hyperbolic as
the generic homogeneous realization — is not what that classification
supports, and is retired.

**Axiom A3 (no privileged directional structure).** The abstract branching
process carries no distinguished structure on its space of directions beyond
a real inner product: its symmetry group is the full orthogonal group
\(O(d)\), not a proper subgroup (not \(U(n)\), not \(Sp(n)\), not a
distinguished orientation). A faithful host is required to realize this
symmetry as basepoint-fixing isometries that are automorphisms of the
Heintze data \((N,A)\).

A3 is an axiom of the same kind as the isotropy switch already in the meter.
It is not implied by the packing bound, Theorem 4.4, or saturation. Those
are silent on *which* Heintze host is realized.

**Theorem 7.1.** Let \(M=N\rtimes\mathbb R\) be a Heintze host (homogeneous,
simply connected, negative curvature). Suppose the point-isotropy
representation, under A3, contains the full irreducible standard
representation of \(O(d)\) on \(\mathfrak n\cong T_oN\), \(d=\dim\mathfrak n
\ge 1\). Then \(N\) is abelian and the derivation \(A\) is a scalar multiple
of the identity: \(M\) is real hyperbolic space \(\mathbb H^{d+1}\).

*Proof sketch.* Left-invariant Heintze data are a nilpotent Lie algebra
\(\mathfrak n\), an inner product, and a derivation \(A\) whose eigenvalues
have positive real part. Basepoint isometries that are automorphisms of this
data act orthogonally on \(\mathfrak n\) and preserve the bracket and \(A\).

A Lie bracket is an \(O(d)\)-equivariant map \(\Lambda^2\mathfrak n\to
\mathfrak n\). Under \(O(d)\), \(\Lambda^2\mathfrak n\) is the irreducible
adjoint representation (or zero if \(d=1\)). By Schur, a nonzero equivariant
map into the standard representation \(\mathfrak n\) can exist only if that
representation appears in \(\Lambda^2\mathfrak n\). Dimensions match only at
\(d=3\). For \(d=1,2\) and \(d\ge 4\) there is no such map, so the bracket
vanishes and \(N\) is abelian. At \(d=3\) the unique \(SO(3)\)-equivariant
candidate is the cross product, which is not \(O(3)\)-equivariant
(orientation-reversing elements change its sign). Thus there is no nonzero
\(O(d)\)-equivariant bracket for any \(d\ge 1\): \(N\) is abelian.

With \(N\) abelian, \(A\in\mathrm{End}(\mathfrak n)\). Equivariance is
commutation with the irreducible real standard representation of \(O(d)\).
Schur's lemma forces \(A\) to be a real scalar. Abelian \(N\) and scalar
\(A\) are Heintze's characterizing data for real hyperbolic space. \(\square\)

The argument excludes complex and quaternionic hyperbolic space because A3
forbids \(U(n)\) and \(Sp(n)\) as the process symmetry. If one assumed only
\(SO(d)\), the \(d=3\) cross-product would survive; full \(O(d)\) is load-
bearing.

**What this closes, and what it does not.** Inside the Heintze class, A3
implies isotropy implies real \(\mathbb H^{d+1}\). That implication is a
theorem. It does not prove that a real generator satisfies A3, that a real
host is homogeneous, or that an embedding cannot add a preferred axis the
process did not have. Those remain honesty items. The former conjecture's
empirical remainder is whether packing rate depends on angle, askable once
growth class is exponential and a sector-wise \(h_{\mathrm{pack}}\) is
defined. Distinguishing exponential from Euclidean occupancy is Corollary
4.3 plus the growth-class gate, not E9 as a forcing test.

---

## 8. The normalized state equation

Raw sectional curvature changes when radial units change. Under
\(d\mapsto a\,d\),

\[
c\mapsto ac,
\qquad
\kappa\mapsto\frac{\kappa}{a^2}.
\]

Therefore the dimensionless curvature per generative step

\[
\bar\kappa:=c^2\kappa
\]

is invariant. The curvature floor becomes

\[
\boxed{
\bar\kappa\ge
\left(
\frac{h_{\mathrm{eff}}\ln2}{n-1}
\right)^2.
}
\]

Its capacity-saturating equality case is the normalized geometric state
equation:

\[
\boxed{
\bar\kappa^*=
\left(
\frac{h_{\mathrm{eff}}\ln2}{n-1}
\right)^2.
}
\]

At \(n=2\),

\[
\boxed{\bar\kappa^*=(h_{\mathrm{eff}}\ln2)^2.}
\]

The familiar formula without \(c\) is therefore a formula for normalized
curvature—or for raw curvature after choosing the process-time gauge \(c=1\).

---

## 9. Beyond isotropy

Without homogeneity and isotropy, there need not be one meaningful sectional
curvature scalar. The coordinate-free bound remains:

\[
\beta\le c\,h_{\mathrm{vol}}.
\]

For pinched negatively curved manifolds, comparison geometry may bracket
\(h_{\mathrm{vol}}\) using curvature bounds. For products, mixed histories, or
strongly anisotropic spaces, volume entropy itself is the state variable; no
single \(\kappa\) should be inferred unless the chosen model justifies it.

The safe general prediction is therefore about total asymptotic capacity, not
about how curvature is distributed direction by direction. A directional
"curvature conservation" law would require additional comparison theorems and
is not assumed here.

---

## 10. Dynamics and stability

The addressability theorem is a static constraint. A non-negative mismatch
function such as

\[
V(\kappa,\kappa^*)
=
(\sqrt\kappa-\sqrt{\kappa^*})^2
\]

is positive-definite and has a unique zero on the physical domain
\(\kappa,\kappa^*\ge0\). It does not,
by itself, prove Lyapunov stability, global attraction, or biological
self-organization.

Those claims require a differentiable potential on an invariant positive
domain and an explicit evolution law, for example

\[
\dot\kappa=-\partial_\kappa U,
\]

plus a proof that \(\dot V<0\) away from equilibrium and an empirical reason to
believe the modeled system follows that dynamics.

Until such a dynamics is supplied, the correct claims are:

- uniqueness of the isotropic equality value;
- monotonic increase of required capacity with rate;
- strict positivity of mismatch away from equality;
- variational selection of the boundary under a stated capacity-cost
  objective.

---

## 11. Measurement protocol

Every test must estimate the two sides independently. The growth-class
gate refuses short radial windows; see [`MEASURABILITY.md`](MEASURABILITY.md).

### Generator side

Estimate \(\beta\) from the generating process:

- substitutions or heritable changes per explicitly defined event;
- sound-change transitions;
- branching events;
- another domain-specific process clock.

The estimator must not use measured curvature.

### Representation side

Estimate \(c\,h_{\mathrm{pack}}\) or \(c\,h_{\mathrm{vol}}\) from the static
metric:

- packing or covering growth;
- tree-ball growth;
- a justified isotropic fit \(c(n-1)\sqrt\kappa\);
- an anisotropic volume-entropy estimator.

The radial calibration \(c\) and its units must be reported.

### Classification side

Measure quartet deviations or Gromov \(\delta\) independently of the state
equation. Treat \(n=2\) as an embeddability prediction for genuinely branching
tree metrics in the selected smooth host class, not as evidence when it has
been back-solved using the equality being tested.

### Report

For each system report

\[
(\beta,\ c,\ h_{\mathrm{vol}},\ \eta,\ \delta,\ n,\ \bar\kappa)
\]

with uncertainties and estimator provenance. Equality should be tested through

\[
\eta=
\frac{\beta}{c\,h_{\mathrm{vol}}},
\]

not produced by defining one axis from the other.

---

## 12. Falsification and scope

### The addressability bound would be challenged by

An operationally demonstrated system that simultaneously:

1. produces retained distinguishable histories at \(\beta>0\);
2. preserves a fixed distinguishability resolution;
3. addresses them at finite radial rate \(c\);
4. inhabits a host with \(c\,h_{\mathrm{pack}}<\beta\).

Under the definitions above, such a counterexample would contradict the
packing count itself.

### Capacity saturation would be challenged by

Systematic, independently measured efficiencies \(\eta<1\) in systems claimed
to optimize description length, after finite-size, unit, distortion, and
anisotropy effects are controlled.

This would not refute the addressability bound. It would refute or restrict the
physical saturation hypothesis.

The ranked experimental protocol — interventions first, then pre-registered
predictions, with named kill lines — is
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md).

### The pure-tree model would be challenged by

Persistent quartet violations, \(\delta>0\), or superior non-tree
representations after matched-capacity comparison. These observations indicate
reticulation or mixing rather than a failure of the packing theorem.

### Outside the theory

The theorem does not constrain systems that do not retain their histories,
systems whose distinguishability resolution collapses, or systems whose
address radius is allowed to grow superlinearly without cost.

---

## 13. Minimal statement

\[
\boxed{
\text{retained exponential novelty}
\Longrightarrow
\beta\le c\,h_{\mathrm{cap}}
}
\]

with \(h_{\mathrm{cap}}=h_{\mathrm{pack}}\) at the fixed operational
resolution. Separately imposing

\[
\boxed{
\beta=c\,h_{\mathrm{cap}}
\quad\text{and}\quad
h_{\mathrm{cap}}=(n-1)\sqrt\kappa
\Longrightarrow
c^2\kappa=
\left(
\frac{h_{\mathrm{eff}}\ln2}{n-1}
\right)^2
}
\]

and independently,

\[
\boxed{
\text{four-point condition}
\Longleftrightarrow
\text{tree metric}
}
\]

For genuinely branching trees in the stated smooth hyperbolic embedding class,
this further supports minimal ambient \(n=2\); path trees remain
one-dimensional. The quartet fixes the qualitative class. The process rate and radial
calibration fix the quantitative scale. The inequality is the enduring limit.
The state equation is its optimal isotropic equality case.
