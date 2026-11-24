import inspect
import sympy.polys.rootisolation as RI
print(inspect.signature(RI.dup_refine_real_root))
exec(open('common.py').read())
from fractions import Fraction as F
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
g11l = [int(c) for c in sp.Poly(G.subs({a1:1,a2:1}), a4).all_coeffs()]
# refine root in (1,2)
lo, hi = RI.dup_refine_real_root(g11l, F(1), F(2), sp.ZZ, steps=40)
print("width 2^-?", (hi-lo), " ~ 2^", sp.log(hi-lo,2).evalf() if hi>lo else 0)
print("lo=", sp.nsimplify(lo), float(lo))
print("hi=", float(hi))
