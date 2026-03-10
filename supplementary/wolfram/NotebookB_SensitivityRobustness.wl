
(* ::Title:: *)
(* Notebook B — Sensitivity & Intervals *)
(* Created: 2025-11-05T07:15:58.907818 *)

$Assumptions = Element[h, Reals] && h>0 && Element[n, Reals] && n>=2;
κ[h_, n_] := (h Log[2]/(n-1))^2;
dκdh = D[κ[h,n], h] // FullSimplify;
dκdn = D[κ[h,n], n] // FullSimplify;
Print["∂κ/∂h = ", dκdh];
Print["∂κ/∂n = ", dκdn];

With[{h0=1.6, σh=0.10, n0=2.0, σn=0.06},
 σκ = Sqrt[(dκdh /. {h->h0,n->n0})^2 σh^2 + (dκdn /. {h->h0,n->n0})^2 σn^2] // N[#, 40]&;
 Print["Propagated σ_κ at (h=1.6±0.1, n=2±0.06): ", σκ];
];

With[{hmin=1.50, hmax=1.70, n0=2.0},
 kLow = κ[hmin, n0] // N[#,50]&; kHigh = κ[hmax, n0] // N[#,50]&;
 Print["κ ∈ [", kLow, ", ", kHigh, "] for h ∈ [1.50,1.70], n=2"];
];
