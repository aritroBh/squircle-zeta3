
import os, json
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
with open('progress.log','a') as fh:
    fh.write("step_02: exact W3(k) k=3..15, 70dps quadrature residual<3e-69 all k; controls k=3 (thm:M3) & k=4 (rem:k4) PASS; "
             "DELTA=-z3/8-lam*pi2/24-lam^3/12 constant; linearity PASS; w3_values.json saved\n")

# Test k=2: predicted affine kernel + predicted W3(2) vs paper thm:M2
import sys
sys.path.insert(0, WORK); sys.path.insert(0, os.path.join(BASE,'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
import sympy as sp
from fastred import *
from buildintegrand import build_integrand_pf
import fastred

pf_int = build_integrand_pf(2)
red = fastred.reduce_k_fast(2, integrand_pf=pf_int)
S2 = sp.sqrt(2); vv = sp.Symbol('v')
def pf_show(P):
    out = {}
    for r, c in P[0].items():
        out[str(r.sympy())] = sp.simplify(c.sympy())
    return out
print("k=2 reduced R_20:", pf_show(red[(2,0)]))
print("k=2 reduced R_11:", pf_show(red[(1,1)]))
# predicted affine formulas at k=2:
pred20 = (1-S2)/(4*vv) + sp.Rational(1,4)/(vv-(1+S2)) + sp.Rational(2*2-7,4)/(vv+(1+S2)) - sp.Rational(2-2,2)/(vv+(S2-1))
pred11 = sp.Rational(2*2-5,4)/vv + S2/(4*vv) - sp.Rational(2*2-5,2)/(vv+(1+S2))
print("R_20(2) matches affine prediction:", sp.simplify(sp.together(pf_to_sympy(red[(2,0)]) - pred20)) == 0)
print("R_11(2) matches affine prediction:", sp.simplify(sp.together(pf_to_sympy(red[(1,1)]) - pred11)) == 0)
