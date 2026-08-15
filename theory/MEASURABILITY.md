# Layer 0 — What a finite point set can decide

## Status

This note opens the measurement layer of the program. Layers I and II are
asymptotic. Every biological sample is finite. The two regimes were being
treated as approximation and limit; they are different objects.

The load-bearing result is a measurability theorem for the only qualitative
question the certified instruments can currently ask: is the host exponential
or polynomial? The theorem takes the finite data a matrix actually supplies —
\(N\) points, radial span ratio \(r\), resolvable shell count \(k\), resolution
\(\varepsilon\) — and returns a three-way verdict:

- **unmeasurable** — no test of the two growth classes has risk \(\le\alpha\);
- **measurable** — a named test is certified at this \((N,r,k,\alpha)\);
- **undecided** — between the bounds.

Unmeasurable datasets are not failed experiments. They are datasets the
growth-class question cannot be asked of. The theorem eliminates those runs
instead of adding a premise to them.

The geometric identities (Lemma 1.1, Proposition 1.2, Proposition 1.3) are
elementary and Lean-checked. The testing bounds (Theorems 2.1 and 2.2) are
paper proofs under an explicit Poisson-increment model. The 2×-span death of
the regression gate is the theorem's empirical shadow, reproduced in
`tests/test_growth_class_gate.py`.

---

## 0. Why this layer exists

Every quantity in Layers I and II is a limsup as radius goes to infinity.
Every measurement is a finite \(\varepsilon\)-separated set. The addressability
bound itself has a finite form — \(N(R)\le P(B(o,r(R)),\varepsilon)\) — and
that finite form is the only quantity in the older spine that lives in the
same regime as the data. Everything else was an asymptotic claim being read
off a finite code.

A theorem of the form "if \(A\) and \(B\) then \(C\)" is a gift to the theory
and a debt to the experiment: to use \(C\) one must first establish \(A\) and
\(B\). Almost every recent theorem in this program lengthened that chain. A
measurability theorem is the rare kind that shortens it. It names the
datasets that cannot answer a question, so those experiments are not run.

This is also the honest finite-blocklength turn. Short codes have different
optimal structure from capacity-achieving ones. Biological representation is
a finite code: every clade is finite in size, finite in depth, and truncated
by loss. "Does biology achieve asymptotic capacity" may be the wrong
question, not because it is hard, but because the asymptotic regime is not
the one the data inhabit. Layer 0 asks what can be known from finitely many
points. Layer I remains the capacity theorem. Layer II remains the
realization theory. Neither is retracted.

---

## 1. Geometric identities

Let \(n(\rho)=|X\cap B(o,\rho)|\) be the occupancy of a finite
\(\varepsilon\)-separated set \(X\) of \(N\) points in a pointed metric space
\((M,d,o)\). Write \(\rho_{\min}\) for the smallest positive observed radius
used in a fit and \(\rho_{\max}\) for the largest, and set

\[
r:=\frac{\rho_{\max}}{\rho_{\min}}>1.
\]

The two model classes compared by the growth-class gate are

\[
n_{\mathrm{exp}}(\rho)=A\,e^{h\rho},
\qquad
n_{\mathrm{poly}}(\rho)=B\,\rho^{d},
\]

with \(h>0\) and \(d>0\). These are composite families. The hardest pair
inside them is the pair that agrees at the two endpoints of the observed
window — the only pair a finite sample is forced to treat as close.

### Lemma 1.1 — endpoint matching

Fix a window \([\rho_{\min},\,r\rho_{\min}]\) and an inner occupancy \(m>0\).
The unique exponential and the unique polynomial that match at both
endpoints satisfy

\[
h=\frac{d}{\rho_{\min}}\cdot\frac{\log r}{r-1},
\qquad
n_{\mathrm{exp}}(r\rho_{\min})
=n_{\mathrm{poly}}(r\rho_{\min})
=m\,r^{d}.
\]

The two rates are not independently choosable once the window is fixed. A
finite sample that only sees span \(r\) cannot tell "rate \(h\)" from
"degree \(d\)" except through the shape of the occupancy *between* the
endpoints.

### Lemma 1.2 — midpoint gap (Lean)

On the scaled interval \(t=\rho/\rho_{\min}\in[1,r]\), the log-occupancy
difference of the endpoint-matched pair is

\[
f(t)
:=
d\log t
-
d\cdot\frac{\log r}{r-1}\,(t-1).
\]

At the geometric midpoint \(t=\sqrt{r}\),

\[
\boxed{
f(\sqrt{r})
=
d\log r\Bigl(\tfrac12-\frac{1}{\sqrt{r}+1}\Bigr).
}
\]

Call \(\gamma(r,d):=\exp\bigl(f(\sqrt{r})\bigr)\) the **midpoint occupancy
ratio**. It depends on \(r\) and \(d\) only. The inner radius, the
resolution, and the sample size have dropped out.

### Proposition 1.3 — exact maximum gap (Lean)

\(f(1)=f(r)=0\), and \(f'(t)=0\) at the single critical point
\(t^\star=(r-1)/\log r\) (the logarithmic mean of \(1\) and \(r\)). The
maximum log-gap is therefore

\[
\boxed{
\Delta(r,d)
:=
d\Bigl(\log\frac{r-1}{\log r}-1+\frac{\log r}{r-1}\Bigr)
=
\max_{t\in[1,r]}|f(t)|.
}
\]

\(\Delta(r,d)/d\) is a function of the span ratio alone. It is the
**span information per dimension**: the most nats of occupancy-shape
information a window of ratio \(r\) can carry between an exponential and a
polynomial of degree \(d\).

### Proposition 1.4 — small-span expansion

As \(r\to 1^+\),

\[
\boxed{
\Delta(r,d)
=
\frac{d}{8}(r-1)^{2}
+
O\bigl((r-1)^{3}\bigr).
}
\]

The signal that distinguishes the two classes vanishes *quadratically* with
the excess span. This is why a modest radial window is not "slightly worse
than a large one." It is a different statistical object.

Values of \(\Delta(r,2)/2\), the span information per dimension:

| \(r\) | 1.5 | 2 | \(e\) | 3 | 4 | 8 |
|---|---|---|---|---|---|---|
| \(\Delta/d\) | 0.0205 | 0.0597 | 0.123 | 0.148 | 0.234 | 0.511 |

At 2× span the entire geometric signal is 0.06 nats per dimension.

---

## 2. The testing problem

### Observation model

The nested occupancies \(n(\rho)\) are a single monotone path; they are not
independent observations. The independent information is in the *increments*
across resolvable shells. Let \(k\) be the number of \(\varepsilon\)-separated
radial shells in the window, so

\[
k
\le
1+\frac{\rho_{\max}-\rho_{\min}}{\varepsilon}
=
1+\frac{\rho_{\min}(r-1)}{\varepsilon}.
\]

**Poisson-increment model.** The occupancy increments on those \(k\) shells
are independent Poisson random variables whose means equal the increment of
the model occupancy. This is an idealization (a deterministic packing is
more regular than Poisson). It is the standard testing model, and extra
regularity would only make distinguishability easier, so an
indistinguishability verdict under Poisson is conservative for the lower
bound's *direction* and is not claimed to be sharp for real packings.

Two simple hypotheses: the endpoint-matched pair of Lemma 1.1, written
\(P_{r,d}\) (polynomial) and \(Q_{r,d}\) (exponential). The composite
classes \(\mathrm{Poly}(d_{\max})\) and \(\mathrm{Exp}(h_{\min})\) contain
this pair, so a lower bound for testing the simple pair is a lower bound
for testing the composite classes.

### Theorem 2.1 — indistinguishability (Le Cam)

Let \(H^{2}(P_{r,d},Q_{r,d})\) be the Hellinger distance of the two
increment laws, equivalently

\[
H^{2}
=
1-\exp\Bigl(-\tfrac12\int_{\rho_{\min}}^{r\rho_{\min}}
\bigl(\sqrt{\lambda_{\mathrm{exp}}}-\sqrt{\lambda_{\mathrm{poly}}}\bigr)^{2}
d\rho\Bigr)
\]

in the continuous-shell limit, where \(\lambda_{\mathrm{exp}}\) and
\(\lambda_{\mathrm{poly}}\) are the endpoint-matched intensities. Then every
(possibly randomized) test of \(P_{r,d}\) against \(Q_{r,d}\) has risk at
least

\[
\boxed{
R^{\star}
\ge
\tfrac12\bigl(1-\sqrt{2H^{2}-H^{4}}\bigr).
}
\]

In particular, if \(R^{\star}\ge\alpha\), the pair is **unmeasurable at
confidence \(\alpha\)**. No improvement of the estimator, and no further
biological modelling, can create a test that is not there.

A convenient closed upper bound on the Hellinger integral, sufficient for
an unmeasurability verdict though not always tight, uses
\(\Delta=\Delta(r,d)\) and the window occupancy \(N_{\mathrm{win}}\):

\[
\sum_{j}(\sqrt{\mu_{j}}-\sqrt{\nu_{j}})^{2}
\le
N_{\mathrm{win}}\bigl(e^{\Delta/2}-1\bigr)^{2}
\]

when the increment log-ratios are themselves bounded by \(\Delta\). The
instrument computes the Hellinger integral of the two intensities directly
and uses this bound only as a check.

### Scaling

Detection requires \(\Delta(r,d)\sqrt{N_{\mathrm{win}}}\) of order one.
Proposition 1.4 then gives the critical excess span

\[
r^{\star}-1
\asymp
N_{\mathrm{win}}^{-1/4}.
\]

The exponent \(-1/4\) is the content. Enlarging the sample is a very weak
way to buy a shorter window. Representative thresholds at \(d=2\),
demanding \(\Delta\sqrt{N}\gtrsim 3\) (a rough composite-test allowance
above the simple-hypothesis Le Cam line):

| span \(r\) | \(\Delta(r,2)\) | \(N\) needed |
|---|---|---|
| 1.5 | 0.041 | \(\sim 5.4\times 10^{3}\) |
| 2 | 0.119 | \(\sim 6.3\times 10^{2}\) |
| 3 | 0.297 | \(\sim 1.0\times 10^{2}\) |
| 4 | 0.468 | \(\sim 41\) |

A few hundred points at 2× span sit on the wrong side of the line. A few
thousand points are required before 1.5× span is even a candidate. This is
the proven form of the empirical finding that the growth-class gate dies at
2× span.

### Theorem 2.2 — composite obstruction at few shells

The regression gate compares two affine models,

\[
\log n = a+b\rho
\qquad\text{versus}\qquad
\log n = a+b\log\rho,
\]

on \(k\) shells, by adjusted \(R^{2}\). Each model has two free parameters.

- If \(k\le 3\), both models interpolate or nearly interpolate. Adjusted
  \(R^{2}\) is undefined or meaningless. The gate cannot run.
- If \(k=4\), each residual has one degree of freedom. The probability that
  the wrong model wins on adjusted \(R^{2}\) remains of order one at the
  \(\Delta(r,d)\) values of a 2× window.
- A sufficient condition for the gate to be a valid level-\(\alpha\) test
  of the composite classes is therefore a *joint* lower bound on
  \((k,\,N_{\mathrm{win}},\,r)\), not a bound on \(N\) alone. The instrument
  refuses unless \(k\ge k_{\min}\) (default 6) **and** Theorem 2.1 does not
  already declare the pair unmeasurable.

Theorem 2.2 is the reason a large leaf set at a single radius is worthless
for growth class: \(r=1\) or \(k=1\), regardless of \(N\). Complete trees
observed only at their leaves are unmeasurable. The occupancy profile needs
interior radii.

### What is not claimed

- The Poisson-increment model is not a theorem about biological sampling.
  Replace it with a more regular packing and the numerical thresholds move;
  the identities of §1 and the \(N^{-1/4}\) scaling do not.
- Theorem 2.1 is a lower bound for the simple endpoint-matched pair. A
  better lower bound for the full composite classes would only enlarge the
  unmeasurable region.
- No finite-sample estimator of the *magnitude* \(h_{\mathrm{pack}}\) is
  claimed. Magnitude estimation is a harder problem (E1 already failed it)
  and is not required to refuse a growth-class call.
- No claim is made that a measurable exponential host is hyperbolic, nor
  that a measurable dataset saturates capacity. Layer 0 decides
  askability, not Layer IIa or IIb.

---

## 3. The measurability predicate

A finite pointed sample \((X,o,\varepsilon)\) is **growth-class measurable
at \((d_{\star},\alpha)\)** when

1. \(k\ge k_{\min}\) resolvable shells lie in the fit window, and
2. the Le Cam risk lower bound of Theorem 2.1 at
   \((N_{\mathrm{win}},r,d_{\star})\) is strictly less than \(\alpha\), and
3. a named test (the adjusted-\(R^{2}\) gate, or the midpoint likelihood
   ratio) has been certified at this grid cell.

Otherwise the sample is **unmeasurable** or **undecided**, and the
instrument must not emit a growth class.

The predicate uses only quantities the distance matrix supplies. It does
not require a process clock, a curvature, a weighted duration, or an
isotropy switch. That is the chain-shortening property: the hypotheses of
the theorem are the measurement itself.

The computational form is `tools/growth_class_gate.py`. The meter
(`tools/addressability_meter.py`) calls it and will not promote an
occupancy slope to a packing entropy, nor place a matrix on the
growth-class axis of the phase diagram, unless the predicate returns
measurable.

---

## 4. Consequence for the empirical program

The growth-class × tree-defect phase diagram remains the honest figure, but
its growth axis is now gated. A matrix that fails Layer 0 is not placed on
that axis. This is the first filter of
[`DECISIVE_EXPERIMENTS.md`](DECISIVE_EXPERIMENTS.md) experiment E0.

It also re-reads the standing IIb failures. Saturation (\(\eta\to 1\)) is
an asymptotic property. A good short code is not a capacity-achieving code.
Finite clades should not be expected to look like the \(R\to\infty\)
equality case, and a failed saturation test on a short-span matrix is not
evidence about capacity filling. This does **not** rescue IIb: it forbids
reading IIb on data that cannot support the question. IIb remains an
asymptotic claim with its own kill lines, to be asked only of samples that
Layer 0 marks measurable and that independently justify a magnitude
estimator.

---

## 5. Formal boundary

| Statement | Status | Artifact |
|---|---|---|
| Lemma 1.2, midpoint identity | THEOREM (Lean) | `theory/lean/ActiveGeometry/Measurability.lean` |
| Proposition 1.3, maximum gap | THEOREM (Lean) | same |
| Proposition 1.4, small-span expansion | THEOREM (paper) | this note |
| Theorem 2.1, Le Cam lower bound | THEOREM (paper; Poisson model) | this note |
| Theorem 2.2, few-shell obstruction | THEOREM (paper; regression) | this note |
| 2×-span death of the \(R^{2}\) gate | EMPIRICAL, reproduced | `tests/test_growth_class_gate.py` |
| Measurability predicate | INSTRUMENT | `tools/growth_class_gate.py` |

Lean does not formalize Hellinger distance, Le Cam's lemma, or the
Poisson-increment model. Those are paper proofs under a named sampling
hypothesis, in the same sense that the relational-capacity lower bound is a
paper proof under Skenderi's theorem.
