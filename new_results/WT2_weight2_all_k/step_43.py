import sys, importlib, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

ok = True
for k in [2,3,4,5,6,7,8]:
    ref = {key: pf_from_sympy(R) for key, R in v_integrand(k).items()}
    mine = build_integrand_pf(k)
    assert set(ref) == set(mine)
    for key in ref:
        d = pf_sub(ref[key], mine[key])
        zero = (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
        if not zero:
            ok = False; print("MISMATCH", k, key)
print("build_integrand_pf == v_integrand for k=2..8:", ok)
