# B2: own T1 descent from MY initial data, level k, down to degree max(6,k-8);  compare (a) T1 vs TOWN1 pickle, (b) candidate(forced) vs T1, (c) candidate vs pickle.
import sys, time, pickle
sys.dont_write_bytecode=True
from core import *
import cand
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
k1, k2 = int(sys.argv[1]), int(sys.argv[2])
def toG(q):
    assert q[0]==0 and q[2]==0, q
    return G(fr(q[1]), fr(q[3]))
t0=time.time()
tot = dict(cells=0, t1_bad=0, cand_bad=0, candpk_bad=0, adj=0, diag=0, off=0, d78=0, cand_adj_k2p=0)
ex=[]
for k in range(k1, k2+1):
    nlow = max(6, k-8)
    T = tower_levelk(k, nlow)
    for key, vals in T.items():
        p,q = key; n=p+q
        if n==k: continue                  # initial level (checked in B1)
        for dl, u in enumerate(vals):
            pk = toG(H[key][dl][k])
            tot['cells']+=1
            if dl>=7: tot['d78']+=1
            if p==q: tot['diag']+=1
            elif p==q+1: tot['adj']+=1
            else: tot['off']+=1
            if not (pk==u):
                tot['t1_bad']+=1
                if len(ex)<8: ex.append(('T1vsPK',key,dl,k))
            c = cand.cand_value(p,q,dl,k)
            if not (c==u):
                tot['cand_bad']+=1
                if len(ex)<8: ex.append(('CANDvsT1',key,dl,k))
            if not (c==pk):
                tot['candpk_bad']+=1
            if p==q+1 and k==2*p: tot['cand_adj_k2p']+=1
    print(f'  level k={k} done  t={time.time()-t0:.0f}s  running={tot}', flush=True)
print('B2', k1, k2, tot, ex, f'time={time.time()-t0:.0f}s')
