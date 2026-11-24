import os
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
result = r"""# LEAD1: the generating function of the leading k-coefficients of the weight parts of M_k

Folder: `research/08_theory/LEAD1_leading_genfun/`.  Everything claimed below was printed by a script in this
folder (`step_*.py`; rerun driver `verify_lead1.py`, ~20 s, ALL PASS).  `mp.mp.dps` set after all imports.
Nothing was written outside this folder.

## VERDICT

**The leading coefficients `c(w) = [k^{w-2}] W_w(k)` for w = 2..10 were computed from the exact top kernel
coefficients `D_{w-1}` by quadrature (405 digits for w <= 9, 250 for w = 10); controls c(2),c(3),c(4),c(6) pass
to ~410 digits.  All of c(5)..c(9) are identified in the Gamma-world monomial basis with integer relations
IDENTICAL at two precisions (dps 420/D=390 and dps 250/D=225).  The generating function is**

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!)  *  [y^w] E(y),
    E(y)  =  Gamma(1 + y/2) Gamma(1/2 - y/2) / sqrt(pi)
          =  2^y Gamma(1-y) Gamma(1+y/2) / Gamma(1-y/2)
          =  exp( y log 2 + sum_{n>=2} s_n y^n / n ),   s_n = zeta(n) (n even),  s_n = (1 - 2^{1-n}) zeta(n) (n odd),

`eps_m = +1,-1,-1,+1` (period 4), `A_n` Euler zigzag numbers (`sec x + tan x = sum A_n x^n/n!`).
Checked against ALL computed c(w): w = 2..9 to 406 digits, w = 10 to 237 digits.

Equivalent forms: with `eps_m A_{m+1}/m! = [z^m](1 - sinh z)/cosh^2 z` (symbolic, step_54),

    F(u) = sum_{w>=2} c(w) u^w  =  (1/2) * [ u^2 (1 - sinh u)/cosh^2 u  \odot  E(u) ]   (Hadamard product),
    sum_{w>=2}  2 (w-2)!/(eps_{w-2} A_{w-1}) * c(w) u^w  =  E(u) - 1 - u log 2.

Proof status (details in genfun.md): the top-coefficient table `D_d` is EXACT for d <= 9 (all keys/poles,
zigzag/binomial pattern; d >= 10 open = tower law T).  GIVEN the table, the reduction of c(w) to the Gamma
factor E is PROVED for even d (odd weights) for all d — parity collapse of the key-sum plus a
generating-function/Beta-integral evaluation — and the single-integral evaluations are PROVED for both
parities for all d.  For odd d (even weights) the identification of the key-sum with the single integral is
checked to 60-80 digits for d <= 9 (real kernel data) but the all-d algebraic collapse is NOT proved
(finite-part obstruction described in genfun.md Sec 3.4).  Hence the law is PROVED given the D-table for odd
weights w (in particular c(3), c(5), c(7), c(9)), and CHECKED at 237-406 digits for all w <= 10.

## 1. Values (60 digits; 405/250 in leading_coeffs.json)

| w | c(w) |
|---|---|
| 2 | 0.531346770191606965284879423243172540237375713200336006150606 |
| 3 | -0.463054519888073265020183366821329463829457953996109311061807 |
| 4 | -0.512151959461527680222200558606808431407661025164465171146113 |
| 5 | 0.40948119421256435691521921506651565742131739799942432899534 |
| 6 | 0.33568755822402035700767227984351182721153530440389179214486 |
| 7 | -0.253122416960440487223947539490993007431297696004986069718015 |
| 8 | -0.189242541413378208209141417583276128090266362460591787429044 |
| 9 | 0.137263340792142938015670190422387507617099134460586704785821 |
| 10 | 0.0984598871418917972706502940864227127442040035844511774009152 |

Sign pattern +,-,-,+,+,-,-,+,+,... (= eps_{w-2}, period 4).

## 2. Controls (step_14, dps 410; two independent quadrature runs, step_27)

| control | residual |
|---|---|
| c(2) - (pi^2 + 6 lam^2)/24 | 0.0 (exact to 410 digits) |
| c(3) + (zeta(3)/8 + pi^2 lam/24 + lam^3/12) | 6.2e-412 |
| c(4) + (pi^4/320 + pi^2 lam^2/48 + lam^4/48 + zeta(3) lam/8) | 1.2e-411 |
| c(6) - (w6 form, 8 monomials/9!) | -2.5e-411 |
| run2 (dps 250) - run1 (dps 420), w=2..9 | <= 6.5e-238 |

## 3. Identifications (Gamma-world monomials; lam = log 2, z_k = zeta(k))

LLL (flint, gram='exact') relations, IDENTICAL at dps 420 (D=390) and dps 250 (D=225); each also EQUALS the
relation derived from the law's closed form (step_25):

    3456 c(5) = 12 lam^5 + 20 pi^2 lam^3 + 180 lam^2 z3 + 9 pi^4 lam + 30 pi^2 z3 + 270 z5
    362880 c(6) = 168 lam^6 + 420 pi^2 lam^4 + 5040 lam^3 z3 + 378 pi^4 lam^2 + 2520 pi^2 lam z3 + 22680 lam z5 + 61 pi^6 + 3780 z3^2
    29030400 c(7) = 1464 lam^7 + 5124 pi^2 lam^5 + 76860 lam^4 z3 + 7686 pi^4 lam^3 + 76860 pi^2 lam^2 z3 + 691740 lam^2 z5
                    + 3721 pi^6 lam + 230580 lam z3^2? -- machine-readable: lll_relations_dps400.json --
    (w = 8, 9: 14 and 19 monomials, max coefficient 2.6e9 and 5.9e10; residuals 5e-397, 3e-395 at dps 400)

Closed forms for all w <= 10 (from the law, matching the LLL relations): `leading_coeffs.json`,
`cw_closed_forms.json`.  Structure: `c(w) = A_{w-1} * (polynomial)/2(w-2)!`, with the top pi-power carrying
extra zigzag factors (61 = A_6 in c(6), 21437 = 17*1261 in c(8), 2449 = 31*79 in c(10)).

## 4. The exact top-coefficient table (D_d), d <= 9

Extracted from `kern_k.pkl` (WTS/theory) by Newton differences over k = d+1..18 (step_08, step_57);
`delta_d = [x^d](sech x + tanh x)`:

    pole -a (all keys):              D_d = (-1)^{j+1} binom(d,j) delta_d
    pole  a (odd d only):            D_d = (-1)^{j}   binom(d,j) delta_d (x 1/2 on the last key)
    pole -b (even d, j < d/2):       D_d = (-1)^{j}   binom(d,j) delta_d
    pole  0 (diagonals):             D_d = (-1)^{j}   binom(d,j) delta_d / 2   (plus -b, odd-d last key)
    everything else:                 0

EXACT for d = 1..8 (all keys, all poles; step_36) and d = 9 (step_57: determined fit on k = 10..18,
zero spare points, still matches the pattern everywhere).  This is the combinatorial heart of the law:
`delta_d` supplies the zigzag factor `eps_{w-2} A_{w-1}` and the endpoint integrals supply `e_w`.

## 5. Evidence for the generating function

| check | where | residual |
|---|---|---|
| law vs c(w), w = 2..9 (dps 420 values) | step_21 | <= 4.8e-406 |
| law rerun vs c(w), w = 2..9 (dps 240) | verify_lead1.py D | <= 9e-228 |
| law vs c(10) from exact D_9 (dps 250) | step_60 | -4.6e-237 |
| E(y) series == Gamma product, y = 0.1 | verify_lead1.py D0 | -1.1e-31 (13-term trunc.) |
| even-d single integral == (d!/2) e_{d+1}, d = 0..16 | step_34/verify E | <= 1e-61 |
| odd-d single integral == (d!/2) e_{d+1}, d = 1..11 | step_41 | <= 2.5e-59 |
| key-sum == single integral, d = 1..8 (real data) | step_50 | <= 5.3e-43 |
| odd-d collapse d = 9 (real D_9 data) | step_61 | 1.6e-37 (dps 80) |
| (F(y)+F(-y))/2 == (E(y)-E(-y))/(2y), y = 0.3 | verify E | -4.3e-46 |
| (piece1 - boundary) - piece2 == (d!/2) e_{d+1}, odd d <= 9 | step_56 | <= 4.2e-46 |

## 6. Files

* `leading_coeffs.json` — values (405/250 digits), closed forms, LLL relations, law statement.
* `genfun.md` — the closed form and its derivation, with the proved/checked status per step.
* `verify_lead1.py` — reruns A (exact D-table), B (quadrature + controls), C (LLL identifications),
  D (law), E (proof-side integral identities); prints PASS/FAIL; ALL PASS in ~20 s.
* `Dtop.json`, `Dtop_d9.json` — exact top coefficients; `Nints_dps420.json`, `Nints_dps250.json` — master integrals;
  `cw_dps420.json`, `cw_dps250.json` — assembled c(w); `lll_relations_dps400.json` — identifications.

## 7. What is NOT proved

* The D_d table for d >= 10 (tower law T of WTS/theory; d = 9 is exact but with zero spare points).
* The odd-d (even-weight) collapse of the key-sum for all d: checked d <= 9, proof blocked by the
  finite-part splitting of G_odd = B(x,y)/2 + sqrt2 2^{-(x+y)} J(x,y) (genfun.md Sec 3.4).
  The even-d case IS proved for all d given the table.
"""
open('result.md','w').write(result)
print('result.md written', len(result))
