(* ::Title:: *)
(* Notebook D — Null Simulations & Falsification (Wolfram) *)
(* Created: 2025-11-05 *)

(* BFS depths computation *)
Clear[BFSDepths];
BFSDepths[adj_List, root_:1] := Module[{n = Length[adj], depths, q, head = 1, tail = 1},
  depths = ConstantArray[-1, n]; q = ConstantArray[0, n];
  depths[[root]] = 0; q[[tail]] = root;
  While[head <= tail,
    With[{u = q[[head]]}, head++;
      Do[
        If[depths[[v]] == -1, depths[[v]] = depths[[u]] + 1; tail++; q[[tail]] = v],
        {v, adj[[u]]}
      ];
    ];
  ];
  depths
];

(* Ball histogram V(R) *)
Clear[BallHist];
BallHist[depths_List] := Module[{Rmax = Max[depths]},
  Table[Count[depths, _?(# <= R &)], {R, 0, Rmax}]
];

(* Fit exponential to estimate κ *)
Clear[FitSlope];
FitSlope[hist_List, minFrac_:0.3] := Module[{Rmax=Length[hist]-1, R, start, x, y, fit, s, ci},
  If[Rmax<5, Return[{Indeterminate, Indeterminate, {Indeterminate, Indeterminate}}]];
  R = Range[0, Rmax]; start = Max[1, Floor[minFrac Rmax]];
  x = N@R[[start+1;;]]; y = N@Log[hist[[start+1;;]] /. 0->1];
  fit = LinearModelFit[Transpose[{x,y}], x, x];
  s = fit["BestFitParameters"][[2]];
  ci = fit["ParameterConfidenceIntervals"][[2]];
  {s^2, s, ci^2}  (* conservative CI by squaring endpoints *)
];

(* Generate b-ary tree *)
Clear[BaryTree];
BaryTree[b_Integer, depth_Integer] := Module[{adj = {{}}, frontier = {1}, nodeId = 2},
  Do[
    frontier = Flatten[Table[
      Do[adj = Append[adj, {}], {b}];
      Do[
        adj[[u]] = Append[adj[[u]], nodeId];
        adj[[nodeId]] = {u};
        nodeId++,
        {b}
      ];
      Range[nodeId - b, nodeId - 1],
      {u, frontier}
    ]];
    frontier = Flatten[frontier],
    {depth}
  ];
  adj
];

Print["Wolfram Language null simulation functions loaded."];
Print["Use BFSDepths, BallHist, FitSlope, and BaryTree for null control tests."];


