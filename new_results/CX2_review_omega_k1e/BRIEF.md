# CX2 (Codex GPT-6): adversarial review of two new results before they go into the paper

Work ONLY in this directory. Python: <project>/.venv/bin/python (sympy, mpmath). Write RESULT.md here.
Try hard to REFUTE each claim; report every gap. Never claim a number you did not print.

## Claim A (ω_n ≠ 0). In the paper (<project>/paper/ieee/main.tex, Proposition "prop:an", search
"\label{prop:an}"), for odd n the coefficient of ζ(n) in a_n is ω_n = (4(−1)^n/n)(1 − 2^{1−n}V_n), V_n = (d/dt)^n √((1+e^{4t})/2) at t = 0,
and the paper says "We do not know whether ω_n is nonzero for every odd n; we checked n ≤ 21". Proposed proof: g(t) = e^t·√(cosh 2t);
z = (cosh 2t − 1)/4 has integer derivatives at 0 and z(0) = 0, cosh 2t = 1 + 4z; √(1+4y) = 1 + Σ_{j≥1} a_j y^j with a_j = (−1)^{j−1}·2·C_{j−1}
(Catalan), all even integers; products of Hurwitz series (integer derivatives at 0) are Hurwitz, so c_k = (√cosh 2t)^{(k)}(0) are integers,
c_0 = 1, c_k even for k ≥ 1; V_n = Σ_k C(n,k) c_k ≡ 1 (mod 2); hence V_n ≠ 2^{n−1} for n ≥ 2 and ω_n ≠ 0.
Check every step (is the composition argument valid for an infinite sum? is the formula for ω_n itself consistent with a_3, a_5 of the
paper? check V_n parity to n = 61 by exact arithmetic) and say whether the proof is complete.

## Claim B (hexagon theorem). ../K1e_hexagon_closed/result.md (read it) claims P'''(1) of every centrally symmetric hexagon in closed form,
with the explicit ζ(3)-coefficient equal to the local Z_e for every ratio r, "canonical" because the trilogarithm functional equations exit
the 11-function basis. Review: (1) is "canonical" correctly stated (as functions of r, are the 11 basis functions linearly independent over
rational functions of r — argue or test numerically at several r); (2) pick a hexagon NOT used by K1e (e.g. vertices (2,0), (3,2), (0,1) and
negatives, your own weights), evaluate K1e's formula (formula_table.txt / Earc_basis_clean.pkl; read step_21–26 for how the parameters are
built) and compare with your own direct computation of P'''(1) from the polar perimeter integral of K_p = {Σ w_j|k_j cross(V_j,x)|^p ≤ 1}
(central finite differences in p with Richardson, ≥ 20 digits). Report agreement or the discrepancy.
