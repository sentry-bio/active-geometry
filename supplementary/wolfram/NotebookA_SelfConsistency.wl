
(* ::Title:: *)
(* Notebook A — Self-Consistency & Uniqueness *)
(* Created: 2025-11-05T07:15:58.907818 *)

$Assumptions = Element[h, Reals] && h>0 && Element[n, Reals] && n>=2 && Element[r, Reals] && r>0 && Element[κ, Reals] && κ>0;

eq = Exp[r Sqrt[κ (n-1)]] == 2^(h r);
sol = Assuming[$Assumptions, Solve[eq, κ, Reals] // FullSimplify];
Print["Closed-form κ(h,n): ", sol];

Assuming[$Assumptions && n==2,
 Module[{F, dF, lim0p, limInf},
  F[kk_] := Sqrt[kk] - h Log[2];
  dF = D[F[κ], κ] // Simplify;
  lim0p = Limit[F[κ], κ->0, Direction->"FromAbove"];
  limInf = Limit[F[κ], κ->Infinity];
  Print["dF/dκ = ", dF];
  Print["F(0+)= ", lim0p, " ; F(∞)= ", limInf];
  Print["Conclusion: strictly increasing with sign change ⇒ exactly one root."];
 ]
];

Assuming[$Assumptions && Element[a, Reals] && a>0 && Element[b, Reals],
 Module[{eqRescaled, kRescaled},
  eqRescaled = Exp[(a r+b) Sqrt[κ (n-1)]] == 2^(h (a r+b));
  kRescaled = Solve[eqRescaled /. r->1, κ, Reals] // FullSimplify;
  Print["Rescaled solution κ'(a,b): ", kRescaled];
 ]
];

h0 = SetPrecision[1.6, 120];
n0 = 2;
kPred = (h0 Log[2]/(n0-1))^2 // N[#, 120]&;
Print["κ( h=1.6, n=2 ) at 120 digits: ", kPred];
