
import os, pickle, sys
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
record = pickle.load(open('record_wt3.pkl','rb'))

# DL Taylor coefficients at v=i: DL = sum_r d_r/(v-r); 1/(v-r) at v=i: 1/(i-r) * 1/(1+w/(i-r)) = sum_n (-1)^n w^n/(i-r)^{n+1}
def taus(DL, n):
    IQ = QSI(0,0,1,0)
    out = ZERO
    for r, d in DL.items():
        out = out + d * ((-ONE)**n) / ((IQ - r)**(n+1))
    return out
t1 = [taus(DL1, n) for n in range(6)]
t2 = [taus(DL2, n) for n in range(6)]
print("tau1:", [str(x) for x in t1[:4]])
print("tau2:", [str(x) for x in t2[:4]])

def xH(L, key, m):
    HL = record[L]['rec'][key]['H'][1]; Hm = record[L-1]['rec'][key]['H'][1]
    a = HL[m-1] if len(HL) >= m else ZERO
    b = Hm[m-1] if len(Hm) >= m else ZERO
    return a - b

def gx(L, key, mp):
    return -xH(L, key, mp+1)/QSI(mp)

# check x_30^(L)[3] recursion for L=5..10
print("\nx_30^(L)[3]: data vs recursion")
for L in range(5, 11):
    lhs = xH(L, (3,0), 3)
    rhs = ZERO
    for mp in range(3, 6):
        rhs = rhs - QSI(3) * gx(L, (4,0), mp) * t1[mp-3]
        rhs = rhs - QSI(1) * gx(L, (3,1), mp) * t2[mp-3]
    print(f" L={L}: data={lhs}  rec={rhs}  match={lhs==rhs}")
