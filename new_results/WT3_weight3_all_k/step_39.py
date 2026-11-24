
import os, pickle, sys
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))
from fastred import *
from fractions import Fraction as Fr

# U(v) = 2v/(v^2+1) - sqrt2/2 as PF: order-1 at +-i coeff 1 (p0=1), const -sqrt2/2
U = [{}, [ONE], [ONE], QSI(0,Fr(-1,2))]

def gm(key, k):
    return rec[(k,key)]['Gm']

for key in [(4,0),(3,1),(5,0),(4,1),(3,2)]:
    d = sum(key)
    print(f"== X_{{{key}}}(d') = Gm^(d') - Gm^(d'-1) = xi * U ?")
    ks = sorted(k for (k,kk) in rec if kk==key and (k-1,kk) in rec)
    for k in ks:
        X = pf_sub(gm(key,k), gm(key,k-1))
        # proportional to U? compare order-1 ipos coeff
        c1 = X[1][0] if X[1] else ZERO
        Xi = pf_scale(U, c1) if not c1.is_zero() else pf_zero()
        ok = pf_sub(X, Xi)
        prop = not (ok[0] or any(not e.is_zero() for e in ok[1]) or any(not e.is_zero() for e in ok[2]) or not ok[3].is_zero())
        print(f"   d'={k}: xi={c1}  {'PROP' if prop else 'not prop: ' + str(X)}")
