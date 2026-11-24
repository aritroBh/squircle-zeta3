# Step 1: exact reduction (fastred.reduce_k_fast) for k = 2..11 from the cached integrands.
# Saves the final pieces (all degrees) to final_<k>.pkl, and prints their exact pole content.
import sys, os, time, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastred
from fastred import *

ks = [int(x) for x in sys.argv[1:]] or list(range(2, 12))
for k in ks:
    out = f"final_{k}.pkl"
    if os.path.exists(out):
        print(f"k={k}: final_{k}.pkl exists, skipping", flush=True); continue
    ipf = pickle.load(open(f"binteg_{k}.pkl", "rb"))
    t0 = time.time()
    res = reduce_k_fast(k, integrand_pf=ipf)
    dt = time.time() - t0
    pickle.dump(res, open(out, "wb"))
    print(f"k={k}: reduce {dt:.2f}s, keys={sorted(res)}", flush=True)
    for key in sorted(res, key=lambda x: (sum(x), x)):
        P = res[key]
        print(f"   key {key} deg {sum(key)}: poles={ [(repr(r), repr(c)) for r,c in P[0].items()] } ipos={len(P[1])} ineg={len(P[2])} const_zero={P[3].is_zero()}", flush=True)
