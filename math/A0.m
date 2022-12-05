(* A0 function from [arxiv:hep-ph/9606211 Eq. (B.5)] *)
(* Arguments are interpreted as squared. *)
A0[m_, _] := 0 /; PossibleZeroQ[m]

A0[m_, q_] := m (Delta + 1 - Log[m/q])
