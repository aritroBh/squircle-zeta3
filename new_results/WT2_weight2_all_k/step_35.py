import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf, sqrt_part_PF
from genintegrand import S_poly

v = sp.Symbol('v', positive=True)
ref = {key: pf_from_sympy(R) for key, R in v_integrand(2).items()}
mine = build_integrand_pf(2)
for key in ref:
    d = pf_sub(ref[key], mine[key])
    zero = (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
    print(key, "match:", zero)
    if not zero:
        # compare numerically to find the factor
        rs = pf_to_sympy(ref[key]); ms = pf_to_sympy(mine[key])
        for xv in (sp.Rational(7,10), sp.Rational(13,10)):
            print("   at", xv, "ref:", sp.N(rs.subs(v,xv),20), " mine:", sp.N(ms.subs(v,xv),20))
