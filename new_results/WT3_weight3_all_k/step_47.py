
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))

# X^L higherpart order-m coeff = H(W^L) - H(W^{L-1}) (for L>=d+1)
def xH(L, key, m):
    d = sum(key)
    if L == d:
        init = record[d]['init'][key]
        ipos = init[1]
        return ipos[m-1] if len(ipos) >= m else ZERO
    HL = record[L]['rec'][key]['H'][1]; Hm = record[L-1]['rec'][key]['H'][1]
    a = HL[m-1] if len(HL) >= m else ZERO
    b = Hm[m-1] if len(Hm) >= m else ZERO
    return a - b

def fit_binom(vals_L, nu):
    """vals_L: dict L->QSI; fit sum_r c_r binom(L, r) for r=0..nu; return [c_r]."""
    Ls = sorted(vals_L)
    L0 = Ls[0]
    ys = [vals_L[L] for L in Ls]
    # finite difference in the step-1 direction on index
    out = []
    cur = ys[:]
    for r in range(nu+1):
        out.append(cur[0])
        cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
    # now out[r] = (Delta^r y)(0) but as function of L0 index; binom(L-L0, r) basis
    return L0, out

for key in [(3,0),(2,1),(4,0),(3,1),(5,0),(4,1),(3,2),(6,0),(5,1),(4,2)]:
    d = sum(key)
    print(f"== key {key} d={d} ==")
    for m in range(2, 2*d-2+1):
        nu = (2*d-4-m)//2
        if nu < 0:
            continue
        vals = {L: xH(L, key, m) for L in range(d+1, 16)}
        L0, coefs = fit_binom(vals, nu)
        topc = coefs[nu]
        # report top coefficient and check the lower coefficients are exact fits
        print(f"  m={m} nu={nu}: top binom-coeff c_{nu} = {topc}   (c's: {coefs})")
