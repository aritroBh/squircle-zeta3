import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
from fastred import *

v = sp.Symbol('v', positive=True)
# Test pf_from_sympy vs sympy apart on k=3 integrands
for key, R in v_integrand(3).items():
    P = pf_from_sympy(R)
    back = sp.simplify(sp.together(pf_to_sympy(P)) - sp.together(R))
    print(key, "roundtrip:", back == 0, " real:", {r.sympy(): c.sympy() for r,c in P[0].items()},
          " ipos:", [e.sympy() for e in P[1]])
