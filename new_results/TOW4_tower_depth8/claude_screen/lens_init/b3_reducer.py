# B3: ACTUAL full-reducer towers (WTS tower_k.pkl, k<=14; no T1 shortcut, built from full integrands binteg_k) vs
#     (a) my rho (initial level), (b) my own T1 descent from my initial data, (c) the forced candidate.
import sys, time, pickle
sys.dont_write_bytecode=True
from core import *
import cand
W='<project>/research/08_theory/WTS_weight_k_structure/theory/'
def toG(q):
    assert q[0]==0 and q[2]==0, q
    return G(fr(q[1]), fr(q[3]))
tot=dict(cells=0, init=0, init_bad=0, t1_bad=0, cand_bad=0, d78=0, adj_k2p=0, nonzero_d78=0)
ex=[]; t0=time.time()
for k in range(7, 15):
    HK = pickle.load(open(W+f'tower_{k}.pkl','rb'))
    T = tower_levelk(k, 6)
    for key in T:
        p,q = key; n=p+q
        for dl in range(0, min(8, 2*n-4)+1):
            r = 2*n-2-dl
            v = HK.get(key, {}).get(r)
            val = toG(v) if v is not None else ZERO
            u = T[key][dl]
            tot['cells']+=1
            if dl>=7:
                tot['d78']+=1
                if not val.iszero(): tot['nonzero_d78']+=1
            if n==k:
                tot['init']+=1
                if not (val==u): tot['init_bad']+=1; ex.append(('init',key,dl,k))
                continue
            if not (val==u):
                tot['t1_bad']+=1; ex.append(('T1',key,dl,k))
            c = cand.cand_value(p,q,dl,k)
            if not (c==val):
                tot['cand_bad']+=1; ex.append(('cand',key,dl,k))
            if p==q+1 and k==2*p: tot['adj_k2p']+=1
    print(f'k={k} done {time.time()-t0:.0f}s {tot}', flush=True)
print('B3', tot, ex[:6])
