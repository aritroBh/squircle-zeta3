import sys, importlib, time, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf

v = None
Rw2pf = [{ZERO: QSI(Fr(1,2)), AQ: QSI(Fr(1,2)), -AQ: QSI(-1), -BQ: QSI(Fr(-1,2))}, [], [], ZERO]

loglines=[]
print(f"{'k':>3} {'build(s)':>9} {'reduce(s)':>10} {'kernel==Rw2':>12}")
for k in range(2, 16):
    fn = f"binteg_{k}.pkl"
    t0=time.time()
    if os.path.exists(fn):
        ipf = pickle.load(open(fn,'rb'))
    else:
        ipf = build_integrand_pf(k)
        pickle.dump(ipf, open(fn,'wb'))
    t1=time.time()
    res = reduce_k_fast(k, integrand_pf=ipf)
    t2=time.time()
    d = pf_sub(res[(1,0)], Rw2pf)
    zero = (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
    msg = f"k={k}: build {t1-t0:.2f}s reduce {t2-t1:.2f}s kernel==Rw2: {zero}"
    print(f"{k:>3} {t1-t0:>9.2f} {t2-t1:>10.2f} {str(zero):>12}", flush=True)
    loglines.append(msg)
    assert zero
with open('progress.log','a') as f:
    f.write("CLOSED-FORM INTEGRAND + FAST REDUCER\n" + "\n".join(loglines) + "\n")
print("ALL k=2..15 EXACT: degree-1 kernel == Rw2")
