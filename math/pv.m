(* A0 function from [arxiv:hep-ph/9606211 Eq. (B.5)] *)
(* Arguments are interpreted as squared. *)
A0[m_, _] := 0 /; PossibleZeroQ[m]

A0[m_, q_] := m (Delta + 1 - Log[m/q])

(* B0 function from [arxiv:hep-ph/9606211 Eq. (B.7)] *)
(* Arguments are interpreted as squared. *)
B0[p_, m1_, m2_, q_] := B0[p, m2, m1, q] /; m1 > m2

B0[p_, m1_, m2_, q_] :=
    Delta + Log[q/m2] /; PossibleZeroQ[p] && PossibleZeroQ[m1 - m2]

B0[p_, m1_, m2_, q_] :=
    Delta + 1 + (m1 Log[q/m1] - m2 Log[q/m2])/(m1 - m2) /; PossibleZeroQ[p]

B0[p_, m1_, m2_, q_] :=
    Module[{eps, fB, s, xp, xm},
           s = p - m2 + m1;
           xp = (s + Sqrt[s^2 - 4 p (m1 - I eps)]) / (2 p);
           xm = (s - Sqrt[s^2 - 4 p (m1 - I eps)]) / (2 p);
           fB[x_] := Log[1 - x] - x Log[1 - 1/x] - 1;
           $Assumptions = { p >= 0, m1 >= 0, m2 >= 0, q >= 0, eps > 0};
           Series[Delta - Log[p/q] - fB[xp] - fB[xm],
                  {eps, 0, 0}] // Normal // Refine
          ];

(* ****** alternative implementation of B0 ****** *)

(* Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104 *)
(* Arguments are interpreted as squared. *)
B02[p_, m1_, m2_, q_] :=
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
