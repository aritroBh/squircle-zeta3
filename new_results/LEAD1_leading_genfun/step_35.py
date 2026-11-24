import json, mpmath, os
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
cw = json.load(open('cw_dps420.json'))
forms = json.load(open('cw_closed_forms.json'))
rel = json.load(open('lll_relations_dps400.json'))
leading = {
 "description": "Leading coefficients c(w) = [k^(w-2)] W_w(k) of the pure-weight parts of the squircle coefficients M_k",
 "values_dps405": cw,
 "closed_forms": forms,
 "lll_relations_dps400": rel,
 "generating_law": "c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * [y^w] E(y), E(y) = Gamma(1+y/2) Gamma(1/2-y/2)/sqrt(pi) = 2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2); eps = +1,-1,-1,+1 (period 4), A = Euler zigzag numbers (sec x + tan x)",
}
json.dump(leading, open('leading_coeffs.json','w'), indent=1)
print('wrote leading_coeffs.json')
result = """# LEAD1: generating function of the leading k-coefficients of the weight parts of M_k

Folder: `research/08_theory/LEAD1_leading_genfun/`.  Everything claimed below was printed by a script in this folder
(`step_*.py`, rerun driver `verify_lead1.py`).  `mp.mp.dps` set after all imports and printed with every residual.

## VERDICT

**The leading coefficients c(w) = [k^(w-2)] W_w(k) for w = 2..9 are computed to 405 digits from the exact top kernel
coefficients D_d by quadrature; all four controls (c(2), c(3), c(4), c(6)) pass to ~410 digits.  All of c(5)..c(9) are
identified in the Gamma-world monomial basis, with relations IDENTICAL at two precisions (dps 420/D=390 and dps 250/D=225).
The generating function is**

    c(w)  =  eps_{w-2} * A_{w-1} / (2 (w-2)!)  *  [y^w] E(y),
    E(y)  =  Gamma(1 + y/2) Gamma(1/2 - y/2) / sqrt(pi)
          =  2^y Gamma(1-y) Gamma(1+y/2) / Gamma(1-y/2)
          =  exp( y log 2 + sum_{n>=2} s_n y^n / n ),   s_n = zeta(n) (n even),  s_n = (1 - 2^{1-n}) zeta(n) (n odd),

with eps_m = +1,-1,-1,+1 (period 4) and A_n the Euler zigzag numbers (sec x + tan x = sum A_n x^n/n!:
1,1,1,2,5,16,61,272,1385,...).  Verified against every c(w), w = 2..9, to 406 digits (step_21).  The prefactor
eps_{w-2} A_{w-1}/(w-2)! equals [z^{w-2}] (1 - sinh z)/cosh^2 z (symbolic, w6/s33; rechecked step_22), so equivalently

    F(u) = sum_{w>=2} c(w) u^w  =  (1/2) [ z^2 (1 - sinh z)/cosh^2 z  \odot  E(z) ] (u)

(Hadamard = coefficientwise product).  The "reduced" generating function is an explicit Gamma product:

    sum_{w>=2} 2 (w-2)! / (eps_{w-2} A_{w-1}) * c(w) u^w  =  E(u) - 1 - u log 2.

Status: the c(w) values and identifications are numerical (two-precision integer relations, residuals below);
the law is CHECKED for w <= 9 (406 digits), not proved for all w.  The exact mechanism is established for d <= 8:
the top kernel coefficient table D_d is exactly the zigzag/binomial table (Section: proof), and GIVEN that table the
reduction of c(w) to the Gamma factor E is a finite exact-linear-algebra + quadrature identity per weight.

## 1. Values (405 digits in leading_coeffs.json; first 50 digits)

| w | c(w) (50 digits) | sign |
|---|---|---|
"""
import mpmath as mp
mp.mp.dps = 80
for w in range(2,10):
    v = mp.mpf(cw[str(w)])
    result += f"| {w} | {mp.nstr(v, 50)} | {'+' if v>0 else '-'} |\n"
result += """
Controls (step_14, dps 410): c(2) - (pi^2+6 lam^2)/24 = 0.0 exactly;  c(3) + S = 6.2e-412;
c(4) - [W_2^2 - lam S - 7 zeta(4)/16] = 1.2e-411;  c(6) - (w6 form) = -2.5e-411.
"""
result += """
## 2. Identifications (Gamma-world monomial basis; lam = log 2, z_k = zeta(k))

Exact relations (LLL at dps 400, D=390; IDENTICAL at dps 250, D=225; step_24, step_26/27; also equal to the
closed-form relations of the law, step_25):

    3456 c(5) = 12 lam^5 + 20 pi^2 lam^3 + 180 lam^2 z3 + 9 pi^4 lam + 30 pi^2 z3 + 270 z5
    362880 c(6) = 168 lam^6 + 420 pi^2 lam^4 + 5040 lam^3 z3 + 378 pi^4 lam^2 + 2520 pi^2 lam z3 + 22680 lam z5 + 61 pi^6 + 3780 z3^2
    29030400 c(7) = 1464 lam^7 + 5124 pi^2 lam^5 + 76860 lam^4 z3 + 7686 pi^4 lam^3 + 76860 pi^2 lam^2 z3 + 691740 lam^2 z5 + 3721 pi^6 lam
                    + 230580 pi^4 lam z3?  -- see lll_relations_dps400.json for machine-readable vectors --
"""
open('result.md','w').write(result)
print('wrote result.md (draft)')
