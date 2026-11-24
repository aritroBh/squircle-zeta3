# B1: my own rho engine (Taylor of S at t*, exact) vs TOWN1 pickle initial data H[key][delta][d], ALL keys d<=52, delta<=min(8,2d-4)
import sys, time, pickle
sys.dont_write_bytecode=True
from core import *
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def toG(q):
    return G(fr(q[1]), fr(q[3])), (q[0]==0 and q[2]==0)
t0=time.time(); cells=0; bad=0; nonpure=0; d78=0; diagcells=0; ex=[]
for (p,q) in sorted(H, key=lambda k:(k[0]+k[1],k)):
    d=p+q
    if d<2: continue
    r = rho_tilde(p,q,8)
    for dl in range(0, min(8,2*d-4)+1):
        v, pure = toG(H[(p,q)][dl][d])
        cells+=1
        if p==q: diagcells+=1
        if dl>=7: d78+=1
        if not pure: nonpure+=1
        if not pure or not (v==r[dl]):
            bad+=1
            if len(ex)<5: ex.append((p,q,dl))
print(f'B1 own-rho vs TOWN1 initial data (level k=d, all keys d<=52): cells={cells} (diag {diagcells}; depth7/8: {d78}) bad={bad} non-sqrt2-pure={nonpure} {ex} time={time.time()-t0:.0f}s')
