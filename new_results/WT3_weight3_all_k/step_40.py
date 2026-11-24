
import os
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
step4 = r'''"""WT3 step 4: instrumented reducer recording (W, HigherPart, Gm) for all keys, k=2..15.
Verifies:
 (W1) degree-2 fixed-point equations for all k:
      k=2: W_20 = R_20^init, W_11 = R_11^init (no feeds)
      k>=3: W_20^(k) = -3 Gm_30 DL1 - Gm_21 DL2,  W_11^(k) = -2 Gm_21 DL1
 (W2) affine Gm law for k=3..15:
      Gm_30^(k) = Gm_30^(3) + (k-3) s30,  Gm_21^(k) = Gm_21^(3) + (k-3) s21,
      s30 = -(1/12)(1 - 2 sqrt2 v/(v^2+1)),  s21 = -3 s30.
Saves record_wt3.pkl.
"""
import sys, os, pickle, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
import fastred
from buildintegrand import build_integrand_pf

def reduce_k_record(k):
    """reduce_k_fast with recording of W (pre-split work item), HigherPart, Gm per key."""
    pf_int = build_integrand_pf(k)
    work = {}
    for key, R in pf_int.items():
        kk, RR = fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    init = dict(work)   # initial integrand (folded keys) at level k
    final = {}
    rec = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if deg == 0:
                final[key] = R
                continue
            simple, higher, const = pf_split(R)
            assert const.is_zero(), f"polynomial part at {key}"
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            rec[key] = {'W': R, 'H': higher}
            if not higher[1] and not higher[2]:
                rec[key]['Gm'] = None
                continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            rec[key]['Gm'] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return final, rec, init

record = {}
t0 = time.time()
for k in range(2, 16):
    final, rec, init = reduce_k_record(k)
    record[k] = {'final': final, 'rec': rec, 'init': init}
    print(f"k={k} recorded ({time.time()-t0:.0f}s)", flush=True)
with open(os.path.join(WORK, 'record_wt3.pkl'), 'wb') as fh:
    pickle.dump(record, fh)

def pf_is_zero(P):
    return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

# ---- W1 verification ----
print("\n== W1: degree-2 fixed-point equations ==")
ok1 = True
for k in range(2, 16):
    rec = record[k]['rec']; init = record[k]['init']
    W20 = rec[(2,0)]['W']; W11 = rec[(1,1)]['W']
    if k == 2:
        a = pf_is_zero(pf_sub(W20, init[(2,0)])); b = pf_is_zero(pf_sub(W11, init[(1,1)]))
        print(f" k=2: W20==init {a}, W11==init {b}")
        ok1 &= a and b
    else:
        Gm30 = rec[(3,0)]['Gm']; Gm21 = rec[(2,1)]['Gm']
        rhs20 = pf_add(pf_scale(pf_mul_DL(Gm30, DL1), Fr(-3)), pf_neg(pf_mul_DL(Gm21, DL2)))
        rhs11 = pf_scale(pf_mul_DL(Gm21, DL1), Fr(-2))
        a = pf_is_zero(pf_sub(W20, rhs20)); b = pf_is_zero(pf_sub(W11, rhs11))
        print(f" k={k}: W20==-3Gm30 DL1 - Gm21 DL2: {a};  W11==-2Gm21 DL1: {b}")
        ok1 &= a and b
print("W1:", "PASS" if ok1 else "FAIL")

# ---- W2 verification ----
print("\n== W2: affine Gm law ==")
s30 = [{}, [QSI(0,Fr(1,12))], [QSI(0,Fr(1,12))], QSI(Fr(-1,12))]
s21 = [{}, [QSI(0,Fr(-1,4))], [QSI(0,Fr(-1,4))], QSI(Fr(1,4))]
G3 = record[3]['rec'][(3,0)]['Gm']; G21 = record[3]['rec'][(2,1)]['Gm']
# cross-check G3, G21 against WT2b record values used in step C1-C4
G3x  = [{}, [QSI(0,Fr(1,8),0,Fr(1,8)), QSI(0,0,0,Fr(1,8)), QSI(0,Fr(-1,12),0,0)],
         [QSI(0,Fr(1,8),0,Fr(-1,8)), QSI(0,0,0,Fr(-1,8)), QSI(0,Fr(-1,12),0,0)], QSI(Fr(1,24),Fr(1,8))]
G21x = [{}, [QSI(0,Fr(-3,8),0,Fr(-1,8)), QSI(0,0,0,Fr(-3,8)), QSI(0,Fr(1,4),0,0)],
         [QSI(0,Fr(-3,8),0,Fr(1,8)), QSI(0,0,0,Fr(3,8)), QSI(0,Fr(1,4),0,0)], QSI(Fr(1,8),Fr(-1,8))]
print(" k=3 Gm30 == WT2b record:", pf_is_zero(pf_sub(G3, G3x)), " Gm21 == WT2b record:", pf_is_zero(pf_sub(G21, G21x)))
ok2 = True
for k in range(3, 16):
    Gm30 = record[k]['rec'][(3,0)]['Gm']; Gm21 = record[k]['rec'][(2,1)]['Gm']
    pred30 = pf_add(G3, pf_scale(s30, Fr(k-3)))
    pred21 = pf_add(G21, pf_scale(s21, Fr(k-3)))
    a = pf_is_zero(pf_sub(Gm30, pred30)); b = pf_is_zero(pf_sub(Gm21, pred21))
    print(f" k={k}: Gm30 affine {a}, Gm21 affine {b}")
    ok2 &= a and b
print("W2:", "PASS" if ok2 else "FAIL")
print("STEP4:", "PASS" if (ok1 and ok2) else "FAIL")
'''
with open('step_04.py','w') as fh:
    fh.write(step4)
print("written")
