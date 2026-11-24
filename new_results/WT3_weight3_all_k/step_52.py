
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))

def rho(d, key, t):
    init = record[d]['init'][key]
    ipos = init[1]
    m = 2*d-2-t
    if m < 1: return ZERO
    return ipos[m-1] if len(ipos) >= m else ZERO

for key in [(3,0),(2,1),(4,0),(3,1),(5,0),(4,1),(6,0),(7,0)]:
    d = sum(key)
    r0 = rho(d, key, 0)
    print(f"key {key} d={d}: beta_t = rho[t]/rho[0]:")
    for t in range(0, 7):
        print(f"   beta_{t} = {rho(d, key, t)/r0}")
