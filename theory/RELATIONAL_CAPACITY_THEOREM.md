# The Relational Capacity Theorem

## Status

This note definitively resolves the problem formerly called Conjecture 4.4.

1. The original **unit-edge, exact-radius** formulation is false as quantified.
2. The natural **weighted generative-clock** formulation is true:

\[
\boxed{
C_{\mathrm{rel}}^{\mathrm{wt}}
(\mathbb H_\kappa^n,c;\varepsilon)
=c(n-1)\sqrt\kappa .
}
\]

The converse is the addressability packing theorem. Achievability follows from
Skenderi's construction of free Bishop--Jones semigroups with critical exponent
arbitrarily close to that of an ambient lattice, together with the
quasi-isometric orbit map proved in the same paper.

The result says something stronger than block achievability:

> In a real hyperbolic host, preserving a branching genealogy costs no
> exponential-order capacity. The relational tax is zero as a supremum,
> without requiring or establishing attainment of the boundary rate.

---

## 1. Why the original conjecture is false

The former Conjecture 4.4 fixed:

- a unit-edge rooted tree \(T\);
- a resolution \(\varepsilon>0\);
- a radial rate \(c>0\);
- \(f(\mathrm{root})=o\);
- the exact finite-depth budget
  \(f(T_R)\subseteq B(o,cR)\) for every integer \(R\);
- \(\varepsilon\)-separation of all distinct vertices.

### Proposition 1.1 -- finite-depth obstruction

If \(0<c<\varepsilon\), every admissible tree is the one-vertex tree.
Consequently its growth is zero, whereas
\(c(n-1)\sqrt\kappa>0\).

### Proof

If the root has a child \(v\), then \(v\in T_1\), so the exact radial budget
gives

\[
d(o,f(v))\le c<\varepsilon.
\]

But the root maps to \(o\), and separation gives

\[
d(o,f(v))
=d(f(\mathrm{root}),f(v))
\ge\varepsilon,
\]

a contradiction. Thus the root has no child. \(\square\)

This is not a failure of hyperbolic capacity. It is a mismatch between a
physical rate \(c\), an arbitrary fixed resolution, and a clock that declares
every edge to consume exactly one unit while forbidding even an additive
startup allowance. The general addressability theorem never required graph
depth to be that clock.

---

## 2. Correct formulation: a weighted generative clock

Let \(T\) be a rooted, locally finite tree. Assign every edge \(e\) a positive
duration \(a(e)>0\). Define

\[
\tau(v):=\sum_{e\subset[\mathrm{root},v]}a(e)
\]

and let \(d_\tau\) be the corresponding weighted path metric. Histories
retained through process time \(R\) are

\[
T_R^\tau:=\{v\in V(T):\tau(v)\le R\}.
\]

A map \(f:V(T)\to\mathbb H_\kappa^n\) is a **weighted relational code at radial
rate \(c\)** and resolution \(\varepsilon\) when:

1. \(f(\mathrm{root})=o\);
2. distinct vertices are \(\varepsilon\)-separated;
3. for some \(D\ge1\) and \(K\ge0\),

   \[
   D^{-1}d_\tau(u,v)-K
   \le d(f(u),f(v))
   \le Dd_\tau(u,v)+K;
   \]

4. for some depth-independent \(A_0\ge0\),

   \[
   d(o,f(v))\le c\,\tau(v)+A_0.
   \]

Define

\[
C_{\mathrm{rel}}^{\mathrm{wt}}
(\kappa,n,c;\varepsilon)
:=
\sup_{(T,a,f)}
\limsup_{R\to\infty}
\frac{\log|T_R^\tau|}{R}.
\]

The weighted clock is not an escape hatch. It is the usual way a process with
non-uniform event costs is converted into physical time. Edge durations are
local, path-additive, and drawn from a finite set in the construction below.

---

## 3. The theorem

### Theorem 3.1 -- relational capacity of real hyperbolic space

For every \(\kappa>0\), \(n\ge2\), \(c>0\), and
\(\varepsilon>0\),

\[
\boxed{
C_{\mathrm{rel}}^{\mathrm{wt}}
(\kappa,n,c;\varepsilon)
=c(n-1)\sqrt\kappa .
}
\]

### External theorem used

Aleksander Skenderi, "Free semigroups of large critical exponent,"
*Journal of Topology*, published 10 July 2026,
DOI [10.1112/topo.70087](https://doi.org/10.1112/topo.70087):

- Theorem 1.1 (Theorem 5.1): free Anosov subsemigroups of a non-elementary
  transverse group have critical exponents arbitrarily close to that of the
  ambient group;
- Theorem 3.1 and Proposition 3.2(4): the Bishop--Jones free-semigroup
  construction and its partition-sum lower bound;
- equation (6.5), in the proof of Theorem 6.1: the orbit map of the semigroup's
  tree is a quasi-isometric embedding into the symmetric space.

Precisely, Skenderi's Theorem 3.1 applies to a discrete non-elementary
convergence group equipped with an expanding coarse cocycle. Theorem 1.1/5.1
applies to a non-elementary \(P_\theta\)-transverse subgroup of a connected
semisimple group with finite center and no compact factors, and to a functional
\(\phi\) proper on distinct sequences:
\(\phi(\kappa(\gamma_j))\to\infty\). A cocompact lattice in
\(\operatorname{Isom}^+(\mathbb H_\kappa^n)\) satisfies these rank-one
hypotheses; the Cartan functional is hyperbolic displacement up to the fixed
metric normalization. Rescaling the standard symmetric-space metric to
sectional curvature \(-\kappa\) rescales critical exponents by
\(\sqrt\kappa\) and preserves quasi-isometry.

### Proof

Write

\[
h:=(n-1)\sqrt\kappa.
\]

#### Upper bound

For a weighted relational code and \(v\in T_R^\tau\),

\[
d(o,f(v))\le cR+A_0.
\]

The image \(f(T_R^\tau)\) is an \(\varepsilon\)-separated subset of
\(B(o,cR+A_0)\). Therefore

\[
|T_R^\tau|
\le P(B(o,cR+A_0),\varepsilon).
\]

Taking logarithms, dividing by \(R\), and using the fixed-resolution packing
entropy of \(\mathbb H_\kappa^n\),

\[
\limsup_{R\to\infty}\frac{\log|T_R^\tau|}{R}
\le ch.
\]

Hence \(C_{\mathrm{rel}}^{\mathrm{wt}}\le ch\).

#### An \(\varepsilon\)-separated ambient lattice orbit

Choose a torsion-free cocompact lattice
\(\Lambda<\operatorname{Isom}(\mathbb H_\kappa^n)\). Such lattices exist for
every \(n\ge2\). Its critical exponent is the ambient volume entropy:

\[
\delta(\Lambda)=h.
\]

The set

\[
F_\varepsilon
:=
\{\gamma\in\Lambda\setminus\{1\}:d(o,\gamma o)<\varepsilon\}
\]

is finite by proper discontinuity. The group \(\Lambda\) is finitely generated
and linear, hence residually finite. There is therefore a finite-index subgroup
\(\Lambda_\varepsilon<\Lambda\) containing none of
\(F_\varepsilon\). For distinct \(\gamma,\eta\in\Lambda_\varepsilon\),

\[
d(\gamma o,\eta o)
=d(o,\gamma^{-1}\eta o)
\ge\varepsilon.
\]

Finite index preserves cocompactness and critical exponent, so

\[
\delta(\Lambda_\varepsilon)=h.
\]

#### Bishop--Jones semigroup

Fix any \(0<a<h\). Apply Skenderi's Theorem 3.1 to
\(\Lambda_\varepsilon\), using the rank-one Busemann cocycle whose magnitude is
\[
\|\gamma\|=d(o,\gamma o).
\]

It produces a finitely generated free semigroup
\(\mathcal T=\langle S\rangle^+\subset\Lambda_\varepsilon\). Its orbit is
\(\varepsilon\)-separated, and equation (6.5) gives a quasi-isometric orbit map
from its regular word tree into \(\mathbb H_\kappa^n\).

Proposition 3.2(4), applied at the identity, gives

\[
\sum_{s\in S}e^{-a\,d(o,so)}\ge1.
\]

Assign every edge labelled \(s\) the duration

\[
a_s:=\frac{d(o,so)}{c}>0.
\]

Only finitely many durations occur. Consequently the weighted and unit word
metrics are bi-Lipschitz equivalent, so Skenderi's quasi-isometric orbit map
remains quasi-isometric for \(d_\tau\).

For a word \(w=s_1\cdots s_k\), the triangle inequality and invariance of the
hyperbolic metric give

\[
\begin{aligned}
d(o,wo)
&\le
\sum_{i=1}^k d(s_1\cdots s_{i-1}o,s_1\cdots s_io)\\
&=
\sum_{i=1}^k d(o,s_io)
=c\,\tau(w).
\end{aligned}
\]

Thus the radial condition holds with \(A_0=0\).

#### Weighted growth

For the free monoid on \(S\),

\[
\sum_{w\in\mathcal T}e^{-q\tau(w)}
=
\sum_{k\ge0}
\left(\sum_{s\in S}e^{-q a_s}\right)^k.
\]

At \(q=ca\),

\[
\sum_{s\in S}e^{-q a_s}
=
\sum_{s\in S}e^{-a\,d(o,so)}
\ge1.
\]

Hence the series diverges and the weighted-tree growth exponent is at least
\(ca\). Explicitly, because \(S\) is finite and
\(\min_{s\in S}a_s>0\), the weighted tree is locally finite and proper. For
its counting function
\[
N_\tau(R):=\#\{w:\tau(w)\le R\},
\]
the standard abscissa-of-convergence identity gives
\[
\inf\left\{q:\sum_w e^{-q\tau(w)}<\infty\right\}
=
\limsup_{R\to\infty}\frac{\log N_\tau(R)}R.
\]
Therefore divergence at \(q=ca\) implies weighted growth at least \(ca\).
This constructs an admissible relational code of rate at least \(ca\). Since
\(a<h\) was arbitrary,

\[
C_{\mathrm{rel}}^{\mathrm{wt}}\ge ch.
\]

Together with the upper bound,

\[
C_{\mathrm{rel}}^{\mathrm{wt}}=ch.
\qquad\square
\]

---

## 4. Radial form

The same result can be stated without choosing a process clock. Let \(T\) range
over rooted, locally finite unit-edge trees with tree metric \(d_T\), and let
\(f:V(T)\to\mathbb H_\kappa^n\) satisfy:

- \(f(\mathrm{root})=o\);
- \(f(V(T))\) is \(\varepsilon\)-separated;
- \(f\) is a quasi-isometric embedding for \(d_T\).

Define

\[
N_f(\rho):=
\#\{v:d(o,f(v))\le\rho\}
\]

and define

\[
h_{\mathrm{rel}}^{\mathrm{rad}}
(\mathbb H_\kappa^n,o;\varepsilon)
:=
\sup_{(T,f)}
\limsup_{\rho\to\infty}
\frac{\log N_f(\rho)}{\rho}.
\]

The same packing upper bound and Bishop--Jones lower bound give

\[
\boxed{
h_{\mathrm{rel}}^{\mathrm{rad}}
(\mathbb H_\kappa^n,o;\varepsilon)
=(n-1)\sqrt\kappa .
}
\]

The radial and process-time capacity values differ by the factor \(c\). The
lower-bound construction realizes that scaling with the local durations
\(a_s=d(o,so)/c\); for a general weighted code, the one-sided radial condition
does not identify clock balls literally with ambient radial balls.

---

## 5. Radial concentration near capacity

The capacity theorem has a converse rigidity statement: a code operating near
capacity cannot hide a positive fraction of its histories deep inside the
available radial ball.

Write \(h=(n-1)\sqrt\kappa\), let

\[
N(R):=|T_R^\tau|,
\]

and, for \(0<\delta<1\), define the radially deficient histories

\[
E_\delta(R)
:=
\left\{
v\in T_R^\tau:
d(o,f(v))\le(1-\delta)c\,\tau(v)+A_0
\right\}.
\]

At fixed resolution there is a constant
\(A=A(\varepsilon,n,\kappa)>0\) such that

\[
P(B(o,\rho),\varepsilon)\le A e^{h\rho}
\qquad(\rho\ge0).
\]

### Theorem 5.1 -- finite-radius concentration bound

For every \(R>0\),

\[
\boxed{
\frac{|E_\delta(R)|}{N(R)}
\le
\frac{A e^{hA_0}e^{h(1-\delta)cR}}{N(R)}.
}
\]

If, for constants \(B\ge1\) and additive capacity deficit
\(\Delta_{\mathrm{cap}}\ge0\),

\[
N(R)\ge B^{-1}e^{(ch-\Delta_{\mathrm{cap}})R},
\]

then

\[
\boxed{
\frac{|E_\delta(R)|}{N(R)}
\le
ABe^{hA_0}
e^{-(h\delta c-\Delta_{\mathrm{cap}})R}.
}
\]

Thus the deficient fraction vanishes exponentially whenever
\(\Delta_{\mathrm{cap}}<h\delta c\). This additive deficit is not the
dimensionless efficiency \(\eta=\beta/(ch)\) used elsewhere in the program;
when the growth rate is \(\beta\),
\(\Delta_{\mathrm{cap}}=ch-\beta=ch(1-\eta)\).

### Proof

Every \(v\in E_\delta(R)\) satisfies \(\tau(v)\le R\), hence

\[
d(o,f(v))
\le(1-\delta)cR+A_0.
\]

The images are \(\varepsilon\)-separated, so

\[
|E_\delta(R)|
\le
P(B(o,(1-\delta)cR+A_0),\varepsilon)
\le
A e^{hA_0}e^{h(1-\delta)cR}.
\]

Divide by \(N(R)\); the second statement follows by its assumed lower bound.
\(\square\)

### Limsup and limit forms

Let

\[
\beta:=\limsup_{R\to\infty}\frac{\log N(R)}R
\ge ch-\Delta_{\mathrm{cap}}.
\]

For every \(\xi>0\), there is an unbounded sequence \(R_j\) along which

\[
N(R_j)\ge e^{(\beta-\xi)R_j}.
\]

Along that sequence,

\[
\frac{|E_\delta(R_j)|}{N(R_j)}
\le
A e^{hA_0}
e^{-(h\delta c-\Delta_{\mathrm{cap}}-\xi)R_j}.
\]

If the growth limit exists (or the corresponding lower exponential bound
holds eventually), the same estimate holds for all sufficiently large \(R\).
No shell-regularity hypothesis is needed: the finite-radius packing estimate
uses the full clock ball in both numerator and denominator.

### Capacity concentration

For a capacity-rate code (\(\Delta_{\mathrm{cap}}=0\)), every fixed
\(\delta>0\) has an exponentially vanishing deficient fraction along a
sequence realizing the limsup. If the growth limit exists (or an eventual
capacity-rate lower bound is assumed), this holds for all sufficiently large
\(R\).

Take the probability measure to be uniform on the finite clock ball
\(T_R^\tau\). Assign any value to the root's undefined radial ratio; its mass
vanishes. Weighted clock balls are finite here because their
\(\varepsilon\)-separated images lie in compact hyperbolic balls. For every
fixed \(L\), \(N(L)/N(R)\to0\) along the applicable sequence. Outside that
vanishing set, \(\tau(v)>L\), and the radial upper bound gives

\[
\frac{d(o,f(v))}{c\tau(v)}
\le1+\frac{A_0}{cL}.
\]

Letting first \(R\to\infty\) and then \(L\to\infty\) yields

\[
\boxed{
\frac{d(o,f(v))}{c\tau(v)}
\longrightarrow 1
\quad\text{in probability under the uniform measure on }T_R^\tau.
}
\]

For a fixed positive deficit \(\Delta_{\mathrm{cap}}\), this stronger
\(1+o(1)\) conclusion does **not** follow. Instead, for every margin \(m>0\),

\[
\Pr_{v\sim\mathrm{Unif}(T_R^\tau)}
\left[
\frac{d(o,f(v))}{c\tau(v)}
\le
1-\frac{\Delta_{\mathrm{cap}}}{hc}-m
\right]
\longrightarrow0
\]

along the applicable sequence, provided
\(0\le\Delta_{\mathrm{cap}}<hc\). In terms of dimensionless efficiency this
threshold is \(1-\Delta_{\mathrm{cap}}/(hc)=\eta\).

### Sharpness construction (sketch)

The exponent arithmetic \(h\delta c-\Delta_{\mathrm{cap}}\) is the one forced
by packing and total growth. A candidate showing sharpness, up to arbitrarily
small exponent error, is the rooted wedge of two Bishop--Jones relational
codes placed in separated boundary sectors:

1. an outer code of growth \(ch-\Delta_{\mathrm{cap}}\), obtained by uniformly slowing a
   capacity-approaching weighted clock;
2. an inner code with radial rate \((1-\delta)c\) and growth
   \(h(1-\delta)c\).

The second component lies in \(E_\delta(R)\), while the first dominates the
total count when \(\Delta_{\mathrm{cap}}<h\delta c\). Their count ratio has
exponent

\[
h(1-\delta)c-(ch-\Delta_{\mathrm{cap}})
=-(h\delta c-\Delta_{\mathrm{cap}}).
\]

To upgrade this arithmetic to a theorem requires a sector-gluing lemma:
sector-localized capacity approximators, uniform cross-component separation,
bounded cross Gromov products (hence a global quasi-isometric lower bound), and
two-sided counting asymptotics. Those ingredients are plausible but are not
proved here. Accordingly this paragraph is a sharpness construction sketch,
not an established corollary.

---

## 6. Exact boundary of the result

### Settled

- The former universally quantified unit-edge conjecture is false.
- Radial relational capacity equals hyperbolic volume entropy.
- Weighted-clock relational capacity equals
  \(c(n-1)\sqrt\kappa\).
- Therefore the relational tax

  \[
  \Gamma
  :=
  c(n-1)\sqrt\kappa-C_{\mathrm{rel}}^{\mathrm{wt}}
  \]

  is **zero at exponential order** in real hyperbolic space.

### Not settled

If every source-tree edge is forced to consume exactly one clock unit, the
exact finite-depth budget \(f(T_R)\subseteq B(o,cR)\) is retained, and no
additive startup allowance is permitted:

- the statement is false for \(c<\varepsilon\);
- the case \(c\ge\varepsilon\) is not implied by the theorem above.

That equal-edge problem is a stronger synchronization constraint, not the
addressability capacity of the host. E2's circle-filling construction probes
this stricter subclass and finds endpoint distortion; it does not challenge
the weighted/radial theorem.

## 7. Formal boundary

The packing upper bound is already machine-checked in
`lean/ActiveGeometry/Packing.lean`. The lower bound depends on Skenderi's deep
Bishop--Jones semigroup theorem, residual finiteness of finitely generated
linear groups, and the critical exponent of cocompact hyperbolic lattices.
These cited results are not currently formalized in Mathlib, so the lower bound
is a paper proof, not a Lean theorem. The finite sub-ball count underlying
Theorem 5.1 is machine-checked in `Packing.lean`; its full exponential
asymptotic corollary remains a paper theorem.
