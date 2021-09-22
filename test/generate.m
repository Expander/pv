Generate3[f_, p2_, m1_, m2_, q2_, prec_] :=
    Module[{val = f[p2, m1, m2, q2]},
           {
               N[#, prec]& @ {p2, m1, m2, q2},
               {Re[N[val, prec]], Im[N[val, prec]]}
           }
    ]

GenerateGridData[f_, {min_, max_, step_},  prec_] :=
    Flatten /@ Join @@ ParallelTable[
        Generate3[f, p2, m1, m2, 1, prec],
        {p2, min, max, step},
        {m1, min, max, step},
        {m2, min, max, step}
    ]

ExportData[f_, prec_] :=
    Module[{data, filename = ToString[f] <> ".txt", step = 1/10},
           Print["Generating data for " <> ToString[f]];
           data = GenerateGridData[f, {0, 4, step}, prec];
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]

precision = 17;
Delta = 0;

ExportData[B0, precision];
