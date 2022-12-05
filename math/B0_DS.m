(* Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104 *)
(* Arguments are interpreted as squared. *)

B0[p_, m1_, m2_, q_] := B0[p, m2, m1, q] /; m1 > m2

B0[p_, m1_, m2_, q_] :=
    Delta + Log[q/m2] /; PossibleZeroQ[p] && PossibleZeroQ[m1 - m2]

B0[p_, m1_, m2_, q_] :=
    Delta + 1 + (m1 Log[q/m1] - m2 Log[q/m2])/(m1 - m2) /; PossibleZeroQ[p]

B0[p_, m1_, m2_, q_] :=
    Module[{a = x / s, b = y / s, delta = a - b},
           (Delta - Log[s / q] + 2
            - (1 + delta)/2 Log[a]
            - (1 - delta)/2 Log[b]
            - 2 Omega[a, b])
    ];

Omega[a_, b_] :=
    Omega[a, b, (a + b)/2 - (a - b)(a - b)/4 - 1/4];

Omega[a_, b_, C_] :=
    Module[{sC = Sqrt[C]},
           sC (ArcTan[(1 + a - b) / (2 sC)] + ArcTan[(1 - a + b) / (2 sC)])
    ] /; C > 0;

Omega[a_, b_, C_] :=
    Module[{sC = Sqrt[-C]},
           sC / 2 Log[((a + b - 1)/2 - sC) / ((a + b - 1)/2 + sC)]
    ] /; C < 0;

Omega[a_, b_, C_] := 0 /; PossibleZeroQ[C];
