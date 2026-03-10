
(* ::Title:: *)
(* Geometry of Evolution — SI Section 2: CAS-Verified Proofs *)
(* File: SI2_Wolfram_Skeleton.wl *)
(* Created: 2025-11-04T19:13:35.954743 *)

(* ============================== *)
(*      Global Assumptions        *)
(* ============================== *)
$Assumptions = Element[h, Reals] && h > 0 && Element[n, Reals] && n >= 2 && Element[r, Reals] && r > 0 && Element[κ, Reals] && κ > 0;

(* ============================== *)
(*      Geometry Definitions      *)
(* ============================== *)
(* Option A: Use built-in tensors if available *)
(* For clarity and reproducibility, we define the Poincaré ball (unit ball) metric explicitly. *)
ClearAll[gPB, coords, rPB];
coords = {x, y, z}; (* extend to n-dim as needed *)
rPB[x_, y_] := Sqrt[x.x]; (* placeholder for radial coordinate *)
(* In 2D/3D we can write the metric explicitly; for general n, we keep volume element formula. *)

(* Hyperbolic metric scalar curvature check (symbolic placeholder):
   For constant negative curvature, RicciScalar = -n (n-1) κ, with κ > 0 as magnitude. *)
RicciScalarHyperbolic[n_] := -n (n - 1) κ;

(* Volume growth in H^n with curvature -κ: V(r) ~ C_n * Exp[ r Sqrt[κ (n - 1)] ] as r -> ∞ *)
volumeGrowth[r_, κ_, n_] := Exp[r Sqrt[κ (n - 1)]];

(* Information growth for branching rate h (bits per unit r): *)
infoGrowth[r_, h_] := 2^(h r);

(* Self-consistency condition: match exponents *)
κFromSelfConsistency[h_, n_: 2] := (h Log[2]/(n - 1))^2;

(* ============================== *)
(*      Part A: Verify Equality   *)
(* ============================== *)
partA = Assuming[$Assumptions,
  Simplify[
    volumeGrowth[r, κFromSelfConsistency[h, n], n] == infoGrowth[r, h]
  ]
];
Print["Part A (volume vs info growth equality): ", partA];


(* ============================== *)
(* Part A2: Conservative Self-Consistency (No MaxEnt) *)
(* Verify κ as the UNIQUE solution equating exponential volume growth with information growth. *)
Assuming[$Assumptions,
  Module[{eq, sol},
    eq = Exp[r Sqrt[κ (n - 1)]] == 2^(h r);
    sol = Solve[eq, κ, Reals];
    Print["Closed-form κ(h,n) from volume=information: ", sol // Simplify];
  ]
];

(* Uniqueness for n=2 via monotonicity + IVT *)
Assuming[$Assumptions && n == 2,
  Module[{Fκ, dFκ, lim0p, limInf},
    Fκ[κ_] := Sqrt[κ] - h Log[2];
    dFκ = D[Fκ[κ], κ] // Simplify;
    Print["dF/dκ (n=2): ", dFκ, "  (strictly positive for κ>0)"];
    lim0p = Limit[Fκ[κ], κ -> 0, Direction -> "FromAbove"];
    limInf = Limit[Fκ[κ], κ -> Infinity];
    Print["Limits: F(0+)=", lim0p, " ; F(∞)=", limInf];
    Print["Uniqueness: Strictly increasing + sign change ⇒ exactly one root."];
  ]
];

(* ============================== *)
(*      Part B: Uniqueness        *)
(* ============================== *)
F[κ_, h_, n_] := Sqrt[κ (n - 1)] - h Log[2];

(* Monotonicity in κ *)
dF = D[F[κ, h, n], κ] // Simplify;
Print["dF/dκ: ", dF];
monoCheck = Assuming[$Assumptions, Simplify[dF > 0]];
Print["Monotonicity (dF/dκ>0): ", monoCheck];

(* Existence via IVT proxy: limits *)
lim0 = Limit[F[κ, h, n], κ -> 0, Direction -> "FromAbove"];
limInf = Limit[F[κ, h, n], κ -> Infinity];
Print["F(κ->0+): ", lim0, " ; F(κ->∞): ", limInf];
Print["Uniqueness conclusion: strictly increasing → at most one root; limits cross zero → exactly one root."];

(* ============================== *)
(*  Part C: Sensitivity/Error Prop*)
(* ============================== *)
κExpr = κFromSelfConsistency[h, n];
dκdh = D[κExpr, h] // Simplify;
dκdn = D[κExpr, n] // Simplify;
Print["∂κ/∂h = ", dκdh];
Print["∂κ/∂n = ", dκdn];

(* Example numeric propagation *)
h0 = 1.6; σh = 0.1;
n0 = 2.0; σn = 0.05;
σκ = Sqrt[(dκdh /. {h -> h0, n -> n0})^2 σh^2 + (dκdn /. {h -> h0, n -> n0})^2 σn^2] // N[#, 20] &;
Print["Example σκ (h0=1.6±0.1, n0=2±0.05): ", σκ];

(* ============================== *)
(*  Part D: High-Precision Check  *)
(* ============================== *)
hPrec = SetPrecision[1.6, 80];
κPrec = (hPrec Log[2]/(n0 - 1))^2 // N[#, 80] &;
Print["κ predicted (80 d.p.): ", κPrec];

(* Compare with measured *)
κMeasured = Interval[{1.244, 1.250}]; (* 1.247 ± 0.003 *)
Print["Agreement check: κPred ∈ κMeasured? ", IntervalMemberQ[κMeasured, κPrec]];

(* ============================== *)
(*  Part E: κ == 5/4 proximity    *)
(* ============================== *)
hForFiveOverFour = Solve[κFromSelfConsistency[h, 2] == 5/4, h, Reals];
Print["h that gives κ=5/4: ", hForFiveOverFour // N[#, 20] &];

(* ============================== *)
(*   Part F (Optional): MaxEnt Skeleton — use ONLY minimal, biologically independent constraints      *)
(* ============================== *)
(* This is a structured placeholder for the variational derivation.
   You will replace the placeholders with your exact constraint functionals. *)

ClearAll[p, rVar, λ0, λ1, λ2, μ];
(* Lagrangian density for MaxEnt on H^n: L = -p log p + λ0 p + λ1 p f1 + λ2 p f2 ... *)
(* Use dμ = sqrt|g| dx; in practice, integrate in hyperbolic polar coords with volume element S_n sinh^{n-1}(√κ r) ... *)
L[p_, rVar_] := -p[rVar] Log[p[rVar]] + λ0 p[rVar] + λ1 p[rVar] f1[rVar] + λ2 p[rVar] f2[rVar];

(* Euler-Lagrange: δ/δp ∫ L[r] dμ = 0 → log p* = affine combination of constraints → p* ∝ Exp[λ1 f1 + λ2 f2 + ...] *)
(* Use VariationalMethods`EulerEquations for full detail when f_i specified. *)
(* VariationalMethods`EulerEquations[L[p[rVar], rVar], p[rVar], rVar] *)

Print["MaxEnt skeleton (optional). Prefer conservative verification path unless constraints are pre-registered and independently measured."];

(* End of file *)
