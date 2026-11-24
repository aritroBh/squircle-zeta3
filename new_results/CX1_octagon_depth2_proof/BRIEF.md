# CX1 (Codex GPT-6): PROVE the two depth-2 reductions of the regular octagon (upgrade the octagon closed form from numerical to proved)

Work ONLY in this directory (you may read anything under research/08_theory/). Use the project
python: <project>/.venv/bin/python (mpmath, sympy, python-flint). No network needed. Never claim a number you
did not print. Write your final report to RESULT.md in this directory (statements, proofs, every check with its residual, rerun commands).

## Context
OCT6 (../OCT6_octagon_closed/result.md, reductions.pkl, step_10/11/19 and o6util.py) found P'''(1) of the regular octagon in closed form.
Two of its inputs are only NUMERICAL (LLL at 745 digits):
- R3: Λ = Li_{2,1}(1−√2, −√2/2) = Σ_{n>m≥1} (1−√2)^n (−√2/2)^m /(n² m)  reduces to single Li₃'s + products (the relation is in reductions.pkl);
- R4: the swap Li_{2,1}(1−√2/2, −√2) reduces likewise.
(The A-atom reductions R1, R2 are already proved by three differentiated functional identities.) Note x2 = 1−√2/2 = x1/(x1−1) with
x1 = 1−√2 (Landen pair), and (1−√2)(−√2/2) = 1 − √2/2 + ... — inspect the products of the arguments: the reductions presumably come from a
depth-2 functional equation specialised at a Landen-paired point.

## Task
1. Print R3 and R4 exactly from reductions.pkl and re-check them at 200 digits.
2. PROVE each: find a proven functional equation that specialises to it, e.g. view Li_{2,1}(x, y) as a function of one variable along a
   curve through the octagon point (with the other argument tied by the Landen relation), differentiate, reduce the derivative to depth 1
   (dilog/log) by known identities, integrate, and fix the constant at a convenient point; or use the stuffle Li₂(x)Li₁(y) =
   Li_{2,1}(x,y) + Li_{1,2}(y,x) + Li₃(xy) plus a shuffle/duality relation (e.g. the depth-2 inversion or the Landen transform of
   Li_{1,1}, Li_{2,1}); any classical depth-2 identity is acceptable if cited. Every intermediate identity must be (a) derived symbolically
   and (b) checked at ≥ 60 digits at several points of its domain.
3. Report for R3 and R4: PROVED (with the proof) or not, and what exactly is still missing.
