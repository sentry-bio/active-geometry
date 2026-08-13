# The Addressability Kernel

This is the minimal mathematical core of Active Geometry. It contains one
general inequality and one conditional equality. Everything else is either a
realization of the host capacity, a classification theorem, or an empirical
claim.

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

This is the kernel's only prohibition. In particular, positive retained growth
at finite radial rate requires positive exponential host capacity.

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

### The open theorem

The kernel has exactly one open mathematical problem, posed precisely as
Conjecture 4.4 of [`MATHEMATICAL_SPINE.md`](MATHEMATICAL_SPINE.md): the
relational tree capacity of \(\mathbb H_\kappa^n\) — the supremum growth rate
of rooted trees admitting depth-uniform quasi-isometric codes at radial rate
\(c\) and fixed resolution — equals \(c(n-1)\sqrt\kappa\). The converse half is
Theorem 4.1. The achievability half is open. Its proof would be a coding
theorem for genealogy; its refutation would make the gap below
\(c(n-1)\sqrt\kappa\) a new host invariant. Neither outcome disturbs the
inequality or the block identity.

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
- that tree thinness measures curvature magnitude.

Tree classification by the four-point condition is compatible with the kernel
but logically independent of it. It can motivate a hyperbolic host class or a
minimal embedding dimension; it cannot supply \(\beta\), \(c\), saturation, or
the curvature scale.

## 9. Formal boundary

`lean/ActiveGeometry/Packing.lean` formalizes the finite-depth metric packing
count using Mathlib's canonical `Metric.packingNumber`. In the convergent-rate
case it proves

```text
faithful_representation_addressable :
  Addressable β c hpack
```

from independent limits for represented-history growth, radial rate, and exact
ball-packing growth. The limsup formulation in the full spine is the more
general paper theorem; the Lean theorem deliberately uses ordinary finite
limits to keep the formal kernel minimal.

The same file now proves finite-block achievability:
`exists_optimal_blockCode` constructs a separated finite codebook whose
cardinality equals the exact packing count, and
`exists_optimal_blockCode_of_properSpace` discharges local finiteness. Together
with `card_le_packingCount`, these declarations machine-check the finite-radius
operational-geometric identity. Lean does not yet formalize its asymptotic
limsup corollary or achievability for nested, causal, or relation-preserving
code families.

The finiteness side condition `HasFinitePacking` is not an extra assumption for
the intended host class: `hasFinitePacking_of_properSpace` proves it for every
proper metric space (ℝⁿ, hyperbolic space, and every complete Riemannian
manifold via Hopf–Rinow). Retention is recorded structurally, and represented
counts are proved nondecreasing in depth.

`lean/ActiveGeometry/Addressability.lean` then formalizes the algebra
downstream. Its central predicates mirror the dependency structure:

- `Addressable β c hcap`;
- `CapacitySaturated β c hcap`;
- `IsotropicHyperbolic hcap n κ`.

The theorem `normalized_state_equation` requires the latter two predicates
explicitly. It does not derive them.
