import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
R = v_integrand(2)[(2,0)]
num, den = sp.fraction(sp.together(R))
print("den factored:", sp.factor(den, extension=S2))
ap = sp.apart(sp.together(R), v, extension=[S2, sp.I])
print("sympy apart:")
for t in sp.Add.make_args(ap): print("   ", sp.simplify(t))
P = pf_from_sympy(R)
print("mine: real:", {r.sympy(): c.sympy() for r,c in P[0].items()})
print("mine: ipos:", [e.sympy() for e in P[1]])
print("diff mine - sympy:", sp.simplify(sp.together(pf_to_sympy(P)) - sp.together(ap)))
