
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
PAPER="<project>/paper/squircle_note/scripts/"
sys.path.insert(0,PAPER); sys.path.insert(0,WTS)
from fractions import Fraction as Fr
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
from ibp_reduce import reduce_k
def pf_is_zero(P): return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()
for k in [3,4]:
    slow = reduce_k(k, verbose=False)
    fast = fastred.reduce_k_fast(k, integrand_pf=build_integrand_pf(k))
    # compare ALL canonical keys, not just degree 2
    allkeys = sorted(set(slow.keys()) | set(fast.keys()))
    bad=0
    for key in allkeys:
        sPF = fastred.pf_from_sympy(slow.get(key,0))
        fPF = fast.get(key, pf_zero())
        if not pf_is_zero(pf_sub(sPF,fPF)): bad+=1; print("  mismatch",k,key)
    print(f"k={k}: slow ibp_reduce == fast reducer, ALL {len(allkeys)} keys, {bad} bad")
