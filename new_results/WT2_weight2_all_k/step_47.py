import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

N = 12
# reduced kernels for all k <= N
red = {}
for k in range(1, N+1):
    fn = f"binteg_{k}.pkl"
    ipf = pickle.load(open(fn,'rb')) if os.path.exists(fn) else build_integrand_pf(k)
    red[k] = reduce_k_fast(k, integrand_pf=ipf)

# assemble series for each key: coeff of eps^k
keys = sorted(set().union(*[set(red[k]) for k in red]))
def qsi_str(c):
    return str(c.sympy())
for key in [(1,0),(2,0),(1,1),(3,0),(2,1)]:
    print(f"==== key {key}")
    for k in range(1, N+1):
        if k < sum(key): continue
        P = red[k].get(key)
        if P is None: continue
        terms = []
        for r, c in sorted(P[0].items(), key=lambda rc: float(rc[0].sympy())):
            terms.append(f"[{c.sympy()}]@({r.sympy()})")
        if P[1]:
            terms.append("IPOLE:" + ",".join(str(e.sympy()) for e in P[1]))
        print(f"   k={k}: " + " ".join(terms))
