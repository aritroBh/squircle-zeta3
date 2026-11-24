"""(a) T1 vs recorded towers at depths 3,4 (and all <=4): u^{(k)}(delta) == Lambda(u^{(k)} parents), from DATA only.
   (b) full closed-form candidate U (closed-form rho, e) vs recorded towers: all keys d<=14, all k<=14, delta<=4 (diag k=2p -> rho_pp closed form)."""
import sys, os, time
ORIG = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, ORIG); sys.path.insert(1, '.')
import towlib as TL
from towlib import QSI, ZERO, keys_of_degree
import t1_direct as H
def toG(q):
    assert q.p1 == 0 and q.p3 == 0, q
    return H.G(q.p0, q.p2)
t0 = time.time()
# (a)
bad = 0; n = 0; byd = {}
for k in range(5, 15):
    for dd in range(3, k):          # parents of degree dd+1 <= k
        for key in keys_of_degree(dd):
            p, q = key
            for delta in range(0, min(4, 2*dd-4)+1):
                U = lambda kk, dp, k=k: TL.u(kk[0], kk[1], k, dp)
                lam = TL.Lambda_ij(p, q, k, delta, U)
                n += 1; byd[delta] = byd.get(delta, 0)+1
                if not (lam == TL.u(p, q, k, delta)):
                    bad += 1
                    if bad < 5: print('T1 data mismatch', key, k, delta)
print(f'(a) T1 recursion vs recorded towers: k=5..14, deg 3..k-1, delta<=4: {n} cells, {bad} bad; per depth {byd} ({round(time.time()-t0,1)}s)')
# (b)
bad = 0; n = 0; byd = {}
for k in range(4, 15):
    for dd in range(4, k+1):
        for key in keys_of_degree(dd):
            p, q = key
            if p == q and k < 2*p: continue
            for delta in range(0, min(4, 2*dd-4)+1):
                ue = H.U(p, q, k, delta)
                got = toG(TL.u(p, q, k, delta) * QSI(0, 2, 0, 0))
                n += 1; byd[delta] = byd.get(delta, 0)+1
                if not (ue == got):
                    bad += 1
                    if bad < 5: print('closed-form U mismatch', key, k, delta)
print(f'(b) closed-form tower law (no data used for rho or e) vs recorded towers: d=4..14,k<=14, delta<=4 : {n} cells, {bad} bad; per depth {byd} ({round(time.time()-t0,1)}s)')
