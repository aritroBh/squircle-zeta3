exec(open('common.py').read())
from fractions import Fraction as F
from sympy.polys.rootisolation import dup_isolate_real_roots, dup_refine_real_root
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
A_ = sp.sympify(open('subres1_A.txt').read()); B_ = sp.sympify(open('subres1_B.txt').read())
g11 = sp.Poly(G.subs({a1:1,a2:1}), a4)
coeffs = [int(c) for c in g11.all_coeffs()]
# root in (1,2): refine
iv = ((F(1),F(2)),1)  # ((lo,hi), multiplicity)?? need the right format
# isolate first
rr = dup_isolate_real_roots(coeffs, sp.ZZ)
print("isolated:", rr)
target = None
for item in rr:
    (lo,hi) = item[0]
    if lo >= 1 and hi <= 2 and lo<hi: target=(F(lo.numerator,lo.denominator),F(hi.numerator,hi.denominator))
print("target interval:", target)
# refine to width 2^-60
from sympy.polys.rootisolation import dup_refine_real_root as drr
# dup_refine_real_root expects the poly and interval; use sympy API
import sympy.polys.rootisolation as RI
helpstr = [n for n in dir(RI) if 'refine' in n]
print(helpstr)
