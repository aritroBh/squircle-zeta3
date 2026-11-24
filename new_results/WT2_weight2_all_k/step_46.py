import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

v = sp.Symbol('v', positive=True)
S2 = sp.sqrt(2); a = 1+S2; b = S2-1

def show(k):
    res = reduce_k_fast(k, integrand_pf=build_integrand_pf(k))
    print(f"==== k={k}")
    for key in sorted(res, key=lambda kk: (-sum(kk), kk)):
        P = res[key]
        terms = []
        for r, c in P[0].items():
            terms.append(f"({c.sympy()})/(v-({r.sympy()}))")
        for j_, e in enumerate(P[1], start=1):
            terms.append(f"({e.sympy()})/(v-i)^{j_}")
        if not terms: terms = ["0"]
        print(f"   deg {sum(key)} {key}: " + " + ".join(terms))

for k in [1,2,3]:
    show(k)
