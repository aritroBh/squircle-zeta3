import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

# Check: do any final reduced kernels have poles at +-i? (paper checks k<=5)
bad = []
for k in range(1, 13):
    fn = f"binteg_{k}.pkl"
    ipf = pickle.load(open(fn,'rb')) if os.path.exists(fn) else build_integrand_pf(k)
    res = reduce_k_fast(k, integrand_pf=ipf)
    for key, P in res.items():
        if sum(key) >= 1 and (P[1] or P[2]):
            bad.append((k, key))
print("reduced kernels with surviving +-i poles (k<=12):", bad if bad else "NONE")

# Also check degree-0 key (0,0): is it always zero?
for k in range(1, 9):
    fn = f"binteg_{k}.pkl"
    ipf = pickle.load(open(fn,'rb')) if os.path.exists(fn) else build_integrand_pf(k)
    res = reduce_k_fast(k, integrand_pf=ipf)
    P = res.get((0,0))
    print(f"k={k}: (0,0) piece:", "absent" if P is None else (P[0], P[1], P[3].sympy()))
