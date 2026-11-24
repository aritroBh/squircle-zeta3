
import os, pickle, sys
BASE = '<project>')
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))
from math import comb, factorial

# initial data: record[d]['init'][key] is the folded initial integrand at level d.
# its +i principal part (order 1..2d-2): rho[m] = init ipos[m-1]
IQ = QSI(0,0,1,0); SQ = QSI(0,1,0,0)
def K(d):
    odd = 1
    for m in range(1, 2*d-2, 2): odd *= m
    return (IQ**(d+1)) * QSI(Fr((-1)**d * odd, 2*factorial(d))) / SQ
def rho_top(p, q):
    d = p+q
    return K(d) * ((-ONE)**p) * QSI(comb(d,p) * (1 + (1 if p>q else 0)))

# check t=0 against data and print t=1,2 for keys to spot pattern
for d in range(3, 10):
    for i in range(d, (d-1)//2, -1):
        j = d - i
        init = record[d]['init'][(i,j)]
        ipos = init[1]
        top = 2*d-2
        r0 = ipos[top-1] if len(ipos) >= top else ZERO
        r1 = ipos[top-2] if len(ipos) >= top-1 else ZERO
        r2 = ipos[top-3] if len(ipos) >= top-2 else ZERO
        ok0 = (r0 == rho_top(i,j))
        print(f"d={d} ({i},{j}): t0 match={ok0}  t1={r1}  t2={r2}")
    print()
