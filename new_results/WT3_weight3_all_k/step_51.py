
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))
IQ = QSI(0,0,1,0)

def gcomb(n, k):
    r = Fr(1)
    for t in range(k):
        r *= Fr(n - t, t + 1)
    return r

def rho(d, key, t):
    init = record[d]['init'][key]
    ipos = init[1]
    m = 2*d-2-t
    if m < 1: return ZERO
    return ipos[m-1] if len(ipos) >= m else ZERO

def betas(d, key, umax):
    r0 = rho(d, key, 0)
    bs = []
    for tt in range(0, umax+1):
        rt = rho(d, key, tt)
        acc = rt/r0
        for uu, b in enumerate(bs):
            acc = acc - b * QSI(gcomb(2-2*d, tt-uu)) * (IQ**(uu-tt))
        bs.append(acc)
    return bs

for key in [(3,0),(2,1),(4,0),(3,1),(5,0),(4,1),(6,0),(5,1)]:
    d = sum(key)
    bs = betas(d, key, 8)
    print(f"key {key} d={d}:")
    for u, b in enumerate(bs):
        print(f"   beta_{u} = {b}")
