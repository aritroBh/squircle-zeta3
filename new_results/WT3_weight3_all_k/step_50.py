
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
from math import comb
record = pickle.load(open('record_wt3.pkl','rb'))
IQ = QSI(0,0,1,0)

# rho_ij^(d)[m] for m = 2d-2-t from init principal parts
def rho(d, key, t):
    init = record[d]['init'][key]
    ipos = init[1]
    m = 2*d-2-t
    if m < 1: return ZERO
    return ipos[m-1] if len(ipos) >= m else ZERO

def beta(d, key, u):
    """beta_u = alpha_u/alpha_0 extracted recursively:
    rho[2d-2-t]/rho_top = sum_{u=0..t} beta_u binom(2-2d, t-u) (2i)^(u-t)."""
    bs = []
    for tt in range(u+1):
        r0 = rho(d, key, 0)
        rt = rho(d, key, tt)
        val = rt/r0
        # subtract u'<tt contributions
        acc = val
        for uu, b in enumerate(bs):
            acc = acc - b * QSI(comb(2-2*d, tt-uu)) * (IQ**(uu-tt))
        # coefficient of beta_tt is binom(2-2d,0)*(2i)^0 = 1
        bs.append(acc)
    return bs[u]

for key in [(3,0),(2,1),(4,0),(3,1),(5,0),(4,1)]:
    d = sum(key)
    print(f"key {key} d={d}: beta_u for u=0..6:")
    for u in range(0, 7):
        print(f"   beta_{u} = {beta(d, key, u)}")
