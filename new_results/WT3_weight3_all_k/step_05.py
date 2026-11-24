"""WT3 step 5: tail-tower level decomposition and full-order recursion (exact, k<=15).
X_ij^(L) := W_ij^(L) - W_ij^(L-1) (level-L contribution), GmX likewise.
(a) verify the FULL fixed-point system (*) at PF level for all keys, k=2..15:
    W_ij = delta_{k,d} R_ij^init - (i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_ii DL1)
(b) verify the source-free X-recursion at PF level for all keys, levels 3..15 (where defined):
    X_ij^(L) = -(i+1) GmX_{i+1,j} DL1 - [i>j](j+1) GmX_{i,j+1} DL2 - [j=i-1] i F(GmX_ii DL1)   (L>d)
(c) polynomial degrees of x_ij^(L)[m] (+i order-m coeff of HigherPart of X) in L.
"""
import sys, os, pickle
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
import fastred

record = pickle.load(open(os.path.join(WORK, 'record_wt3.pkl'), 'rb'))

def pf_is_zero(P):
    return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

def feed(key, rec_k):
    """the feed part - (i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_ii DL1)
       from a dict of Gm's (rec level k)."""
    i, j = key
    out = pf_zero()
    GmA = rec_k.get((i+1, j), {}).get('Gm')
    if GmA is not None:
        out = pf_add(out, pf_scale(pf_mul_DL(GmA, DL1), Fr(-(i+1))))
    if i > j:
        GmB = rec_k.get((i, j+1), {}).get('Gm')
        if GmB is not None:
            out = pf_add(out, pf_scale(pf_mul_DL(GmB, DL2), Fr(-(j+1))))
    if j == i-1:
        GmD = rec_k.get((i, i), {}).get('Gm')
        if GmD is not None:
            out = pf_add(out, pf_scale(pf_fold(pf_mul_DL(GmD, DL1)), Fr(-i)))
    return out

# (a) full system (*)
bad = []
for k in range(2, 16):
    rec = record[k]['rec']; init = record[k]['init']
    for key in rec:
        d = sum(key)
        W = rec[key]['W']
        rhs = feed(key, rec)
        if d == k:
            rhs = pf_add(rhs, init[key])
        if not pf_is_zero(pf_sub(W, rhs)):
            bad.append((k, key))
print("(a) full fixed-point system (*) at PF level, k=2..15:", "ALL MATCH" if not bad else bad[:8], flush=True)

# (b) source-free X recursion
bad2 = []
X = {}     # X[(L,key)] = {'W','H','Gm'} differences
GmX = {}
for L in range(3, 16):
    for key in record[L]['rec']:
        d = sum(key)
        if d >= L:   # X at top level = initial data; not source-free
            continue
        WL = record[L]['rec'][key]['W']; Wm = record[L-1]['rec'][key]['W']
        XW = pf_sub(WL, Wm)
        HL = record[L]['rec'][key]['H']; Hm = record[L-1]['rec'][key]['H']
        XH = pf_sub(HL, Hm)
        GL = record[L]['rec'][key]['Gm']; GmL_ = record[L-1]['rec'][key]['Gm']
        XG = pf_sub(GL, GmL_) if (GL is not None and GmL_ is not None) else None
        X[(L,key)] = {'W': XW, 'H': XH, 'Gm': XG}
        # source-free recursion: XW == feed from GmX at degree d+1
        recX = {kk: {'Gm': (X.get((L,kk)) or {}).get('Gm')} for kk in record[L]['rec'] if sum(kk)==d+1}
        rhs = feed(key, recX)
        if not pf_is_zero(pf_sub(XW, rhs)):
            bad2.append((L, key))
print("(b) source-free X-recursion at PF level, L=3..15:", "ALL MATCH" if not bad2 else bad2[:8], flush=True)

# (c) polynomial degrees of x_ij^(L)[m] in L  (HigherPart +i coeffs of X)
print("(c) x-polynomial degrees nu_x(d,m) [claim: floor((2d-4-m)/2) for m<=2d-4; 0 at source]:")
for key in [(3,0),(2,1),(4,0),(3,1),(2,2),(5,0),(4,1),(3,2)]:
    d = sum(key)
    per = {}
    for L in range(d+1, 16):
        if (L,key) not in X: continue
        for m, c in enumerate(X[(L,key)]['H'][1], start=1):
            if m < 2: continue
            per.setdefault(m, {})[L] = c
    line = []
    for m in sorted(per):
        # max poly degree over 4 components
        deg = -1
        for ci in range(4):
            vals = [per[m][L] for L in sorted(per[m])]
            comp = [getattr(c, f'p{ci}') for c in vals]
            cur = comp[:]; r = -1 if all(x==0 for x in comp) else 0
            while len(cur) > 1 and any(x != 0 for x in cur) and r < 12:
                cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]; r += 1
            deg = max(deg, r)
        claim = max(-1, (2*d-4-m)//2)
        line.append(f"m{m}:deg{deg}(cl{claim}){'!' if deg!=claim else ''}")
    print(f"  {key} (d={d}): " + " ".join(line), flush=True)
