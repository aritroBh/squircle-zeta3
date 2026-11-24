import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
from fastred import *

v = sp.Symbol('v', positive=True)
key, R = list(v_integrand(3).items())[1]
print("R =", R)
num, den = sp.fraction(sp.together(R))
print("num =", sp.expand(num))
print("den =", sp.factor(den, extension=S2))
P = pf_from_sympy(R)
diff = sp.together(pf_to_sympy(P)) - sp.together(R)
print("diff numeric:", sp.N(diff.subs(v, sp.Rational(7,10)), 20), sp.N(diff.subs(v, sp.Rational(19,10)), 20))
# compare against sympy apart
ap = sp.apart(sp.together(R), v, extension=[S2, sp.I])
print("apart =")
for t in sp.Add.make_args(ap): print("   ", sp.simplify(t))
