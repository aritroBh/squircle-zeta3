"""Step 1 (exact, no floating point): run the fast reducer for each k and save the final pieces.

Usage: nice -n 10 python w4_reduce.py 2 3 4 5 6 7 8
Writes red_<k>.pkl = dict {(i,j): PF}  (all degrees i+j = 0..k) and prints the structure of the pieces.
"""
import sys, os, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fastred
from fastred import *

ks = [int(x) for x in sys.argv[1:]] or [2, 3, 4, 5, 6, 7, 8]
for k in ks:
    out = os.path.join(HERE, f"red_{k}.pkl")
    if os.path.exists(out):
        print(f"k={k}: red_{k}.pkl exists, skipping", flush=True); continue
    ipf = pickle.load(open(os.path.join(HERE, f"binteg_{k}.pkl"), "rb"))
    t0 = time.time()
    res = reduce_k_fast(k, integrand_pf=ipf)
    t1 = time.time()
    pickle.dump(res, open(out, "wb"))
    print(f"k={k}: reduce {t1-t0:.1f}s", flush=True)
    for key in sorted(res, key=lambda x: (-sum(x), -x[0])):
        P = res[key]
        print(f"   key {key} deg {sum(key)}: real poles {sorted(str(r) for r in P[0])}  ipos {len(P[1])} ineg {len(P[2])} const_zero {P[3].is_zero()}")
