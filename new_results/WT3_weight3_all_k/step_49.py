
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))

def xH(L, key, m):
    HL = record[L]['rec'][key]['H'][1]; Hm = record[L-1]['rec'][key]['H'][1]
    a = HL[m-1] if len(HL) >= m else ZERO
    b = Hm[m-1] if len(Hm) >= m else ZERO
    return a - b

for key in [(2,2),(3,3),(4,4)]:
    d = sum(key)
    print(f"== diagonal {key} d={d}: x^(L)[m] values")
    for m in range(2, 2*d-1):
        vals = {L: xH(L, key, m) for L in range(d+1, 16)}
        if all(c.is_zero() for c in vals.values()): continue
        print(f"  m={m}: ", {L: str(c) for L, c in sorted(vals.items())})
