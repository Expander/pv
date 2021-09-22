Generate3[f_, p2_, m1_, m2_, q2_, prec_] :=
    Module[{val = f[p2, m1, m2, q2]},
           N[#, prec]& @ ComplexExpand @ { p2, m1, m2, q2, Re[val], Im[val] }
    ]

GenerateGridData[f_, {min_, max_, step_},  prec_] :=
    ParallelMap[Generate3[f, Sequence @@ #, 1, prec]&,
                Tuples[Table[x, {x, min, max, step}], 3]]

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
