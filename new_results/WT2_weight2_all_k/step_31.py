import sys, importlib, time, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
S2 = sp.sqrt(2); a = 1+S2; b = S2-1
Rw2 = 1/(2*v) + 1/(2*(v-a)) - 1/(v+a) - 1/(2*(v+b))
Rw2pf = pf_from_sympy(Rw2)

def get_integrand_pf(k):
    fn = f"integ_{k}.pkl"
    if os.path.exists(fn):
        return pickle.load(open(fn,'rb'))
    raw = v_integrand(k)
    d = {key: pf_from_sympy(R) for key, R in raw.items()}
    pickle.dump(d, open(fn,'wb'))
    return d

loglines=[]
for k in [4,5,6,7,8]:
    t0=time.time(); ipf = get_integrand_pf(k); t1=time.time()
    fast = reduce_k_fast(k, integrand_pf=ipf); t2=time.time()
    ker = fast[(1,0)]
    diff = pf_sub(ker, Rw2pf)
    is_zero = (not diff[0]) and (not diff[1]) and (not diff[2]) and diff[3].is_zero()
    msg = f"k={k}: integrand {t1-t0:.1f}s, reduce {t2-t1:.2f}s, kernel==Rw2: {is_zero}"
    print(msg, flush=True); loglines.append(msg)
    assert is_zero
with open('progress.log','a') as f:
    f.write("FAST REDUCER (verified vs original at k=2,3 on all keys)\n" + "\n".join(loglines) + "\n")
print("all k=4..8: kernel == Rw2 EXACT")
