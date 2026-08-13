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

## 3. Saturation is an additional condition

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

## 4. Isotropic hyperbolic realization

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

## 5. Unit-invariant equality

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

## 6. Logical dependency

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

## 7. What is not in the kernel

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

## 8. Formal boundary

`lean/ActiveGeometry/Addressability.lean` formalizes the algebra after
\(\beta\le c\,h_{\mathrm{cap}}\) is supplied as a hypothesis. Its central
predicates mirror the dependency structure:

- `Addressable β c hcap`;
- `CapacitySaturated β c hcap`;
- `IsotropicHyperbolic hcap n κ`.

The theorem `normalized_state_equation` requires the latter two predicates
explicitly. It does not derive them.
