# The Addressability Kernel

This is the minimal mathematical core of Active Geometry. It is organized in
**two layers**, and the whole program's parsimony depends on not confusing
them.

- **Layer I — the universal capacity theory (curvature-free).** Three
  independent quantities, one inequality that always holds, one exact identity
  for the best possible block code, and a ladder of constrained capacities
  beneath it. Nothing here mentions curvature, trees, or biology. It is where
  the generality and the honesty live.

- **Layer II — the curvature realization (where real systems live).** This
  layer has two sublayers that must not be merged, because the evidence for
  them is very different.

  - **Layer IIa — the host class (realization).** *Which* geometry hosts the
    data: is it hyperbolic/tree-like, and at what dimension? This is the
    **better-supported** biological claim — seed-stable embeddings, curvature
    as a fixed design parameter, tree-defect near zero, \(n=2\) as an
    embeddability floor. Occupancy of exponential room is a measurement of a
    process, not a bake-off that forces the host. E9 is a finite-sample
    shadow of Corollary 4.3 at matched packing.

  - **Layer IIb — saturation.** *Whether* a given process fills its budget,
    giving the state-equation equality \(\bar\kappa^*=(h_{\mathrm{eff}}\ln2)^2\).
    This is the **harder, less-supported** claim; every independent test to date
    sits below its own kill line. It requires the extra hypothesis \(\eta=1\).

The one inequality is Layer I. The host class is IIa. The equality is IIb, and
only after both saturation and isotropy are supplied. The four-point tree
theorem is an independent classifier attached to IIa: it can motivate a
hyperbolic host class or a minimal embedding dimension, but it does not
calibrate capacity or curvature.

A warning that the sublayer split makes sharp: **isotropy is asserted, never
measured.** In the meter it is the command-line flag
`assume_isotropic_hyperbolic` — the entire IIa premise reduced to a switch. No
result licenses that switch. A3 and occupancy are the remaining IIa
questions; E9 does not grant isotropy.

---

# Layer I — The universal capacity theory

## 1. Independent quantities

Let \(R\) denote generative depth and let \(N(R)\) be the number of retained,
operationally distinguishable histories through depth \(R\). Define

\[
\beta:=\limsup_{R\to\infty}\frac{\log N(R)}{R}
\quad
\text{[nats / generative step]}.
\]

Let \((M,d,o)\) be the representation space. At a fixed resolution
\(\varepsilon>0\), let

\[
P(\rho,\varepsilon)
:=
P(B(o,\rho),\varepsilon)
\]

be the maximum number of \(\varepsilon\)-separated points in the radius-\(\rho\)
ball. Define the host-capacity rate

\[
h_{\mathrm{cap}}
:=
\limsup_{\rho\to\infty}
\frac{\log P(\rho,\varepsilon)}{\rho}
\quad
\text{[nats / radial distance]}.
\]

Finally, if the histories through depth \(R\) are represented inside a ball of
radius \(r(R)\), define

\[
c:=\limsup_{R\to\infty}\frac{r(R)}{R}
\quad
\text{[radial distance / generative step]}.
\]

The three quantities \(\beta\), \(c\), and \(h_{\mathrm{cap}}\) are logically
independent. No one of them is defined by either of the others.

## 2. The unconditional limit

### Addressability theorem

Assume:

1. distinct retained histories remain at least \(\varepsilon\) apart;
2. their images through depth \(R\) lie in \(B(o,r(R))\);
3. \(0<c<\infty\);
4. \(h_{\mathrm{cap}}<\infty\);
5. the radii enter the asymptotic regime, \(r(R)\to\infty\).

Then

\[
\boxed{\beta\le c\,h_{\mathrm{cap}}.}
\]

### Proof

Faithfulness gives the finite-depth counting inequality

\[
N(R)\le P(r(R),\varepsilon).
\]

Taking logarithms, dividing by \(R\), and applying the definitions of the two
asymptotic rates yields

\[
\limsup_{R\to\infty}\frac{\log N(R)}R
\le
\left(
\limsup_{R\to\infty}\frac{r(R)}R
\right)
\left(
\limsup_{\rho\to\infty}
\frac{\log P(\rho,\varepsilon)}{\rho}
\right).
\]

Hence \(\beta\le c\,h_{\mathrm{cap}}\).

This is the kernel's only prohibition. Written with the quantifiers flipped
it is a trichotomy, not only a ceiling:

- \(\beta>0\) and \(c<\infty\) \(\Rightarrow\) \(h_{\mathrm{cap}}>0\)
  (exponential packing);
- \(h_{\mathrm{cap}}=0\) at fixed dimension \(\Rightarrow\) either
  \(\beta=0\) (forgetting) or \(c=\infty\) (address radius superlinear,
  and exponential in Euclidean degree \(n\));
- every finite-\(R\) census is a packing count on a compact ball; the
  displayed rates are the \(R\to\infty\) limit.

Finite rate is a condition of the first reading, not a standing hypothesis
on every representation. Processes that overwrite rather than retain sit
outside the bound.

## 3. The general block-capacity identity

The converse has a fully general achievability partner when the operational
task is stated at the same level of structure. Define the finite-block
addressing number

\[
A_{\mathrm{block}}(\rho,\varepsilon)
:=
\max\left\{
|C|:\ C\subseteq B(o,\rho),\
C\text{ is }\varepsilon\text{-separated}
\right\}.
\]

In a proper metric space the maximum exists, and by the definition of metric
packing,

\[
\boxed{
A_{\mathrm{block}}(\rho,\varepsilon)
=P(\rho,\varepsilon).
}
\]

Thus the operational capacity to assign distinguishable addresses is exactly
the geometric packing capacity, for every finite radius and every proper
metric host. At a linear radial budget \(\rho=cR\), \(c>0\), monotonicity of
packing counts gives the asymptotic identity

\[
\boxed{
C_{\mathrm{block}}(c,\varepsilon)
:=
\limsup_{R\to\infty}
\frac{\log A_{\mathrm{block}}(cR,\varepsilon)}{R}
=c\,h_{\mathrm{cap}}.
}
\]

No tree, hyperbolicity, stationarity, isotropy, or curvature is needed. This is
the exact operational-geometric unification available in full generality.

Its scope is also exact: a block code may be redesigned at each depth and is
required only to distinguish its terminal messages. It need not form one
nested retained codebook, obey a parent-child motion bound, or preserve a
source hierarchy's pairwise relational metric. Those stronger tasks define
constrained capacities no larger than \(C_{\mathrm{block}}\). Their
achievability requires extension, causal, or embedding hypotheses.

The identity concerns the host's *best possible* code. It does not say that a
given physical hierarchy uses that code or fills it.

### The constrained-capacity ladder

The deepest object of Layer I is not a single number but a **ladder of
capacities**, one per admissibility class \(\mathcal A\) — the meanings of
"represent a hierarchy," ordered by how much structure they demand:

1. **block** — distinguish the terminal messages;
2. **persistent** — retain one nested family of addresses;
3. **causal** — bound parent-to-child motion;
4. **relational** — preserve the source's pairwise metric up to fixed
   distortion (genealogy);
5. reticulate and other structural constraints.

Writing \(C_{\mathcal A}(M,c;\varepsilon)\) for the supremum retained-growth
rate achievable in host \(M\) under class \(\mathcal A\),

\[
\boxed{
C_{\mathrm{relational}}\le C_{\mathrm{causal}}\le C_{\mathrm{persistent}}
\le C_{\mathrm{block}}=c\,h_{\mathrm{cap}}.
}
\]

The top rung is always a theorem (the block identity above). The middle rungs
are, at present, definitions. The relational rung depends on the host: it can
be strictly below block capacity in general, but Theorem 4.4 proves equality
for real hyperbolic space under a local weighted process clock. Active Geometry
is therefore a **constrained capacity theory** — geometry supplies addresses;
the admissibility class decides which addresses are usable — with a complete
relational coding theorem in its canonical host.

### The relational capacity theorem

The problem formerly called Conjecture 4.4 is settled, after correcting a
false clocking condition.

The old exact unit-edge statement was false: if \(c<\varepsilon\), the first
child must be both within \(c\) and at least \(\varepsilon\) from the root.
With a local weighted clock — positive edge durations, path-additive process
time, and a radial budget \(d(o,f(v))\le c\,\tau(v)+A_0\) — the full identity
holds:

\[
\boxed{
C_{\mathrm{rel}}^{\mathrm{wt}}
(\mathbb H_\kappa^n,c;\varepsilon)
=c(n-1)\sqrt\kappa.
}
\]

The upper bound is the packing theorem. The lower bound uses Skenderi's 2026
Bishop--Jones theorem: free semigroup orbit trees have critical exponents
arbitrarily close to the ambient hyperbolic exponent and are
quasi-isometrically embedded. Weighting each generator by its geometric
displacement divided by \(c\) gives the process clock and the radial bound.

Thus genealogy has **zero exponential-order tax** in real hyperbolic space;
the proof does not require or establish endpoint attainment. The proof and the exact boundary of the
stronger equal-edge problem are in
[`RELATIONAL_CAPACITY_THEOREM.md`](RELATIONAL_CAPACITY_THEOREM.md).

The theorems do not decide whether nature instantiates the premises or fills
the budget. The ranked experimental protocol is
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md).

## 4. Saturation of a given hierarchy is an additional condition

Define

\[
\eta:=\frac{\beta}{c\,h_{\mathrm{cap}}}
\]

when the denominator is positive. The theorem gives \(\eta\le1\).

Capacity saturation is the separate condition

\[
\boxed{\beta=c\,h_{\mathrm{cap}}}
\qquad(\eta=1).
\]

It follows, for example, if the represented histories fill the available
address space at exponential order:

\[
\log P(r(R),\varepsilon)-\log N(R)=o(R),
\]

together with convergence of the three relevant rates. It can also be selected
by an optimization problem that minimizes attainable host capacity subject to
lossless addressability. Neither mechanism follows from the inequality alone.

### The slack decomposes

When a process is represented under a relational class, the block slack splits
into two physically distinct terms. With block capacity
\(B:=c\,h_{\mathrm{cap}}\) and relational capacity \(C_{\mathrm{rel}}\),

\[
\boxed{
B-\beta
=
\underbrace{(B-C_{\mathrm{rel}})}_{\text{relational tax }\Gamma}
+
\underbrace{(C_{\mathrm{rel}}-\beta)}_{\text{utilization slack }\Delta_{\mathrm{use}}},}
\qquad
\eta_{\mathrm{block}}
=
\underbrace{\frac{C_{\mathrm{rel}}}{B}}_{\text{availability}}\cdot
\underbrace{\frac{\beta}{C_{\mathrm{rel}}}}_{\text{utilization}} .
\]

The distinction is load-bearing in general hosts. In the weighted hyperbolic
class, Theorem 4.4 proves \(C_{\rm rel}^{\rm wt}=B\), so
\(\Gamma=0\): relation preservation does not reduce exponential capacity.
There the block ratio is the utilization ratio. Outside that class or under
the stronger equal-edge synchronization constraint, an observed
\(\eta_{\rm block}<1\) may still mix relational tax with under-utilization.

### Near capacity forces radial concentration

In \(\mathbb H_\kappa^n\), the same packing count applied to the smaller ball
of histories satisfying
\(d(o,f(v))\le(1-\delta)c\tau(v)+A_0\) gives

\[
\frac{|E_\delta(R)|}{|T_R^\tau|}
\le
C e^{-(h\delta c-\Delta_{\rm cap})R},
\qquad
h=(n-1)\sqrt\kappa,
\]

for an eventual lower growth bound \(ch-\Delta_{\rm cap}\); a limsup gives the
same conclusion along a realizing sequence with an arbitrary exponent margin.
Here \(\Delta_{\rm cap}=ch-\beta=ch(1-\eta)\) is additive deficit, not
dimensionless efficiency. At exact capacity,
\(d(o,f(v))/(c\tau(v))\to1\) in probability under uniform counting on clock
balls, along that sequence (or eventually if the growth limit exists). With
fixed positive deficit, only a one-sided band follows. This is Theorem 5.2 of
the spine: the clock--radius identity is a **conditional consequence of
IIa+IIb**, not an unconditional biological law.

---

# Layer II — The curvature realization

Two sublayers, different evidence. **IIa (host class):** §5 realization and
Theorem 7.1 — *which* geometry hosts the data (hyperbolic? what dimension?),
conditional on axiom A3; the better-supported biological claim. **IIb
(saturation):** the *claim* that a process fills its budget — the
coordinate-free condition \(\eta=1\) is defined in §4, and combined with the
IIa realization it yields the §6 equality; the harder, less-supported claim.
The state equation needs both sublayers.

## 5. Isotropic hyperbolic realization

The capacity rate need not be curvature. Add the independent modeling
assumption that the host is the constant-curvature space
\(\mathbb H_\kappa^n\), with sectional curvature \(-\kappa\),
\(\kappa\ge0\). Its volume entropy is

\[
h_{\mathrm{cap}}=h_{\mathrm{vol}}=(n-1)\sqrt\kappa.
\]

The addressability theorem then gives the curvature floor

\[
\boxed{
\kappa\ge
\left(\frac{\beta}{c(n-1)}\right)^2.
}
\]

If capacity is also saturated, and only then,

\[
\boxed{
\kappa^*=
\left(\frac{\beta}{c(n-1)}\right)^2.
}
\]

Writing \(\beta=h_{\mathrm{eff}}\ln2\) merely converts bits to nats.

### Why curvature, and not just "some exponential host"

Layer I is deliberately curvature-free. Given homogeneity *and isotropy*,
the space-form classification plus polynomial exclusion already select
\(\mathbb H_\kappa^n\); that half is not a conjecture. Among homogeneous
negatively curved hosts, constant-curvature hyperbolic space is the
symmetric special case, not the generic one (Heintze). Theorem 7.1 of the
spine closes the remainder *inside that class*: axiom A3 (full \(O(d)\)
directional symmetry, realized as Heintze automorphisms) forces abelian
\(N\) and scalar contraction, hence real \(\mathbb H^{d+1}\). A3 is
asserted about the generator, not measured; the meter still treats
isotropy as a switch (`--assume-isotropic-hyperbolic`).

Real systems are neither exactly homogeneous nor isotropic, so hyperbolic
is the *natural* chart for an A3 generator, not a forced one.
Reticulation sits near the tree class, with four-point defect \(\delta\)
measuring the distance.

## 6. Unit-invariant equality

Under a radial rescaling \(d\mapsto a d\),

\[
c\mapsto ac,
\qquad
\kappa\mapsto\kappa/a^2.
\]

Therefore

\[
\bar\kappa:=c^2\kappa
\]

is invariant. The conditional equality becomes

\[
\boxed{
\bar\kappa^*
=
\left(
\frac{h_{\mathrm{eff}}\ln2}{n-1}
\right)^2.
}
\]

At \(n=2\),

\[
\boxed{\bar\kappa^*=(h_{\mathrm{eff}}\ln2)^2.}
\]

The unnormalized formula is valid only after the radial gauge \(c=1\) has been
chosen.

## 7. Logical dependency

\[
\begin{array}{c}
\text{faithful finite-rate representation}
\\[2mm]
\Downarrow
\\[2mm]
\beta\le c\,h_{\mathrm{cap}}
\\[4mm]
\begin{array}{cc}
\text{capacity saturation}
&
\text{isotropic hyperbolic realization}
\end{array}
\\[2mm]
\Downarrow
\\[2mm]
\displaystyle
\bar\kappa^*
=
\left(
\frac{h_{\mathrm{eff}}\ln2}{n-1}
\right)^2
\end{array}
\]

The lower statement is not equivalent to the upper statement. It requires both
additional branches.

## 8. What is not in the kernel

The kernel does not prove:

- that a given process retains histories at positive rate;
- that Shannon alphabet entropy estimates \(\beta\);
- that a host is homogeneous, isotropic, or constant-curvature;
- that a physical system saturates capacity;
- that an ambient dimension is selected by optimization;
- that a non-negative mismatch function generates attracting dynamics;
- that tree thinness measures curvature magnitude;
- **axiom A3** (the generator's directional symmetry is full \(O(d)\),
  realized as Heintze automorphisms) — asserted, not measured; Theorem 7.1
  is conditional on it.

Axiom A3 is empirical, not a remaining geometric conjecture. The former
relational-capacity conjecture is now Theorem 4.4. Theorem 7.1 is a paper
sketch conditional on A3, not a Lean theorem. The stronger exact unit-edge
synchronization problem remains unresolved for \(c\ge\varepsilon\), but it is
no longer the definition of host capacity.

The full status ledger — every claim tagged THEOREM / IDENTITY / OPEN /
CONVENTION / INSTRUMENT / EMPIRICAL and bound to its artifact — is
[`CLAIMS.md`](CLAIMS.md), enforced by
[`../tools/check_doc_artifacts.py`](../tools/check_doc_artifacts.py). The
whole-program map, including the two manuscripts and their seams, is
[`PROGRAM.md`](PROGRAM.md).

Tree classification by the four-point condition is compatible with the kernel
but logically independent of it. It can motivate a hyperbolic host class or a
minimal embedding dimension; it cannot supply \(\beta\), \(c\), saturation, or
the curvature scale.

## 9. Formal boundary

`lean/ActiveGeometry/Packing.lean` formalizes the finite-depth metric packing
count using Mathlib's canonical `Metric.packingNumber`. In the convergent-rate
case it proves

```text
addressability_limit :
  Addressable β c hpack
```

from an explicit premise that the radii tend to infinity and independent
ordinary limits (`Tendsto`) for represented-history growth, radial rate, and
exact ball-packing growth, with packing-number finiteness discharged in every
proper metric host. The limsup formulation in the full spine is the more
general paper theorem; the Lean theorem deliberately uses ordinary finite
limits to keep the formal kernel minimal. Theorem 4.4
(Skenderi) and Theorem 7.1 (Heintze / A3) are not this theorem and are not Lean.

The same file proves finite-block achievability:
`exists_optimal_blockCode` constructs a separated finite codebook whose
cardinality equals the exact packing count, and
`exists_optimal_blockCode_of_properSpace` discharges local finiteness. Together
with `card_le_packingCount`, these declarations machine-check the finite-radius
operational-geometric identity. Lean does not formalize the asymptotic
limsup corollary or achievability for nested, causal, or relation-preserving
code families.

The finiteness side condition `HasFinitePacking` is not an extra assumption for
the intended host class: `hasFinitePacking_of_properSpace` proves it for every
proper metric space (ℝⁿ, hyperbolic space, and every complete Riemannian
manifold via Hopf–Rinow). Faithfulness (separation in a growing ball) is the
input to the bound; retention (`points_monotone`) is a separate structure, and
represented counts of a retained representation are proved nondecreasing in
depth.

`lean/ActiveGeometry/Capacity.lean` then formalizes the algebra downstream.
Its central predicates mirror the dependency structure:

- `Addressable β c hcap`;
- `CapacitySaturated β c hcap`;
- `hcap_eq_spaceForm hcap n κ` (rate identification with \((n-1)\sqrt\kappa\),
  not a theorem that the host is \(\mathbb H^n_\kappa\)).

The theorem `normalized_state_equation` in `StateEquation.lean` requires the
latter two predicates explicitly. It does not derive them.
