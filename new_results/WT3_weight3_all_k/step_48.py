
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
from math import comb, factorial
record = pickle.load(open('record_wt3.pkl','rb'))

def xH(L, key, m):
    HL = record[L]['rec'][key]['H'][1]; Hm = record[L-1]['rec'][key]['H'][1]
    a = HL[m-1] if len(HL) >= m else ZERO
    b = Hm[m-1] if len(Hm) >= m else ZERO
    return a - b

def lead_coeff(key, m, nu):
    d = sum(key)
    vals = [xH(L, key, m) for L in range(d+1, 16)]
    cur = vals[:]
    for r in range(nu):
        cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
    # leading binom coefficient = Delta^nu (constant); but should verify constancy
    return cur  # list, should be constant if degree exactly nu

print("z(d,2) leading binom-coefficients (c_{d-3}), key (d-j,j):")
for d in range(3, 9):
    nu = d-3
    row = []
    for j in range(0, (d)//2+1):
        i = d-j
        if i < j: continue
        cur = lead_coeff((i,j), 2, nu)
        row.append(((i,j), cur[0], all(c==cur[0] for c in cur)))
    print(f" d={d}: " + "  ".join(f"{key}:{c}({ok})" for key,c,ok in row))
print("\nz(d,4) leading (c_{d-4}):")
for d in range(4, 9):
    nu = d-4
    row = []
    for j in range(0, (d)//2+1):
        i = d-j
        if i < j: continue
        cur = lead_coeff((i,j), 4, nu)
        row.append(((i,j), cur[0], all(c==cur[0] for c in cur)))
    print(f" d={d}: " + "  ".join(f"{key}:{c}({ok})" for key,c,ok in row))
print("\nz(d,6) leading (c_{d-5}):")
for d in range(5, 9):
    nu = d-5
    row = []
    for j in range(0, (d)//2+1):
        i = d-j
        if i < j: continue
        cur = lead_coeff((i,j), 6, nu)
        row.append(((i,j), cur[0], all(c==cur[0] for c in cur)))
    print(f" d={d}: " + "  ".join(f"{key}:{c}({ok})" for key,c,ok in row))
