
import os, sys, json
from fractions import Fraction as F
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
sys.dont_write_bytecode = True
if WD not in sys.path: sys.path.insert(0, WD)
import qsqrt2
from w4_symbol import Q
import wt6s_core as C6
from wt6s_core import add, symbol_of_kernel_terms, pairing_vector6, ban_float, self_check
ban_float(); self_check()
import sympy as sp
k = sp.Symbol('k')
def parse_poly(s): return sp.expand(sp.sympify(s)) if s is not None else sp.Integer(0)
dd = json.load(open('deg5_residue_polynomials.json'))
# W6 symbol from deg5 polys evaluated at k=5
om = {}
for key in ["50","41","32"]:
    i,j = int(key[0]), int(key[1])
    terms = []
    for pole in ["0","a","-a","-b"]:
        cell = dd[key].get(pole, {})
        rat = parse_poly(cell.get('rat')).subs(k,5); sq = parse_poly(cell.get('sqrt2')).subs(k,5)
        if rat or sq: terms.append([pole,[str(sp.Rational(rat)), str(sp.Rational(sq))]])
    om[f"{i},{j}"] = terms
S = {}
for ij in [(5,0),(4,1),(3,2)]:
    S = add(S, symbol_of_kernel_terms(om[f"{ij[0]},{ij[1]}"], *ij))
pv = pairing_vector6(S)
print("W6 pairing from deg5 POLYNOMIALS extrapolated to k=5:", [str(v) for v in pv])
print("W6 pairing from kern_5.pkl (ACTUAL k=5)              : ['449/4*sqrt2', '-30', '15*sqrt2', '-9']")
# difference
diff = [pv[0]-Q(0,F(449,4)), pv[1]-Q(-30), pv[2]-Q(0,15), pv[3]-Q(-9)]
print("difference (poly - actual):", [str(v) for v in diff])
