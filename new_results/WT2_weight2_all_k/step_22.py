import sys, time, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred
importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
ok = True
for k in [2,3,4]:
    for key, R in v_integrand(k).items():
        P = pf_from_sympy(R)
        d1 = sp.N((sp.together(pf_to_sympy(P)) - sp.together(R)).subs(v, sp.Rational(7,10)), 25)
        d2 = sp.N((sp.together(pf_to_sympy(P)) - sp.together(R)).subs(v, sp.Rational(19,10)), 25)
        good = abs(complex(d1)) < 1e-20 and abs(complex(d2)) < 1e-20
        ok = ok and good
        if not good: print("FAIL", k, key, d1, d2)
print("all conversions OK:", ok)
