
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))
IQ = QSI(0,0,1,0)

def rho(d, key, t):
    init = record[d]['init'][key]
    ipos = init[1]
    m = 2*d-2-t
    if m < 1: return ZERO
    return ipos[m-1] if len(ipos) >= m else ZERO

# t1 = -i(d-2) * rho_top for all keys, d<=15
bad = []
for d in range(2, 16):
    for i in range(d, (d-1)//2, -1):
        j = d-i
        r0 = rho(d, (i,j), 0); r1 = rho(d, (i,j), 1)
        if r1 != r0 * QSI(0,0,-(d-2),0):
            bad.append((d, (i,j), r1, r0 * QSI(0,0,-(d-2),0)))
print("t1 = -i(d-2) rho_top for all keys d=2..15:", "ALL MATCH" if not bad else bad[:6])

# beta_2 for key (d,0): fit closed form
print("\nbeta_2 (key (d,0)): real = -(2d^2-9d+11)/4 ?  imag = ?")
for d in range(3, 13):
    b2 = rho(d, (d,0), 2)/rho(d, (d,0), 0)
    re_pred = QSI(Fr(-(2*d*d-9*d+11), 4))
    print(f" d={d}: beta_2={b2}   re_match={b2.p0==re_pred.p0}")
