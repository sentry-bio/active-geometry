
(* ::Title:: *)
(* Notebook C — Cross-Domain Prediction (Wolfram) *)
(* Created: 2025-11-05T07:15:58.907818 *)
data = {{"Zika",1.54,1.20},{"SARS-CoV-2",1.67,1.32},{"HIV-1",1.77,1.45},{"Measles",1.93,1.58},{"CMV",1.94,1.60},{"All Life",1.60,1.247}};
df = AssociationThread[{"system","h","k_meas"}, #]& /@ data;
kPred[h_] := (h Log[2])^2;
df2 = Map[Append[#, "k_pred"->kPred[#["h"]]]&, df];
r = Correlation[Lookup[df2,"k_pred"], Lookup[df2,"k_meas"]];
Print["Pearson r = ", N[r,20]];
Print["Table: ", df2];
