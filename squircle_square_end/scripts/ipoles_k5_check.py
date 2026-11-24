"""ipoles_k5_check.py -- re-run (not independent) of the exact reducer ibp_reduce.reduce_k at k = 2..5.

Paper: "Investigating zeta(3)", Section IX-F, Step 3 ("for k <= 5 an exact check shows that no pole at
+-i survives the reduction") and Proposition 1. For every reduced piece it tests whether v^2 + 1
divides the denominator, i.e. whether a pole at v = +-i is left; only the real letters 0, +-a, +-b
may remain.
Run from this folder:  python ipoles_k5_check.py 2 3 4 5      (about 2 minutes; reference output in
../outputs/ipoles_k5_check.out)
"""
import os
import sys
import time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v

for k in [int(x) for x in sys.argv[1:]] or [2, 3, 4, 5]:
    t0 = time.time()
    res = reduce_k(k, verbose=False)
    bad = []
    for key, R in res.items():
        num, den = sp.fraction(sp.together(R))
        P = sp.Poly(den, v)
        q, r = sp.div(P, sp.Poly(v**2 + 1, v))      # does v^2 + 1 divide the denominator?
        if r.is_zero:
            bad.append(key)
        res_i = sp.simplify(sp.residue(R, v, sp.I)) if r.is_zero else 0
        print(f"k={k} key={key} weight={sum(key)+1} (v^2+1)|den: {r.is_zero} residue at i: {res_i}", flush=True)
    print(f"k={k}: pieces with a pole at +-i: {bad}   ({time.time()-t0:.0f}s)", flush=True)
