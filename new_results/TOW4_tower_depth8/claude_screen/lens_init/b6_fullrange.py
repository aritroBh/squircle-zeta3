# B6: candidate (forced family from MY a_delta; no T1, no recursion) vs TOWN1 pickle for EVERY level k0<=k<=52, every key d<=52, every delta<=min(8,2d-4).
#     (TOW4's own check_full78 stopped at k<=min(d+8,52).)  k0=d off-diagonal, 2p+1 diagonal; diagonal initial level separately in B1.
import sys, time, pickle
sys.dont_write_bytecode=True
from core import *
import cand
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def toG(q):
    assert q[0]==0 and q[2]==0
    return G(fr(q[1]), fr(q[3]))
t0=time.time(); cells=0; bad=0; d78=0; far=0; ex=[]
for (p,q) in sorted(H, key=lambda k:(k[0]+k[1],k)):
    d=p+q
    if d<4: continue
    for dl in range(0, min(8,2*d-4)+1):
        k0 = d if p>q else 2*p+1
        for k in range(k0, 53):
            if k not in H[(p,q)][dl]: continue
            pk = toG(H[(p,q)][dl][k])
            c = cand.cand_value(p,q,dl,k)
            cells+=1
            if dl>=7: d78+=1
            if k>d+8: far+=1
            if not (c==pk):
                bad+=1
                if len(ex)<6: ex.append((p,q,dl,k))
    if d%10==0 and p==d: print(f'  d={d} t={time.time()-t0:.0f}s cells={cells} bad={bad}', flush=True)
print(f'B6 candidate vs pickle, ALL k up to 52: cells={cells} (depth 7/8: {d78}; beyond TOW4 range k>d+8: {far}) bad={bad} {ex} time={time.time()-t0:.0f}s')
