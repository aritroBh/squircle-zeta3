import sys, importlib, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
ok = True
for k in [2,3,4,5]:
    t0=time.time()
    for key, R in v_integrand(k).items():
        P = pf_from_sympy(R)
        d1 = (sp.together(pf_to_sympy(P)) - sp.together(R)).subs(v, sp.Rational(7,10))
        d2 = (sp.together(pf_to_sympy(P)) - sp.together(R)).subs(v, sp.Rational(19,10))
        good = abs(complex(sp.N(d1,25))) < 1e-20 and abs(complex(sp.N(d2,25))) < 1e-20
        ok = ok and good
        if not good: print("FAIL", k, key, sp.N(d1,20))
print("all conversions OK:", ok)
