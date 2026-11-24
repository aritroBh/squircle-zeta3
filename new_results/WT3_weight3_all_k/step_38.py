
import os, sys
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK); sys.path.insert(0, WORK)
from fractions import Fraction as Fr
from fastred import *
import fastred

# Gm_30^(3), Gm_21^(3) from WT2b record (verified against my reducer below)
G3  = [{}, [QSI(0,Fr(1,8),0,Fr(1,8)), QSI(0,0,0,Fr(1,8)), QSI(0,Fr(-1,12),0,0)],
         [QSI(0,Fr(1,8),0,Fr(-1,8)), QSI(0,0,0,Fr(-1,8)), QSI(0,Fr(-1,12),0,0)], QSI(Fr(1,24),Fr(1,8))]
G21 = [{}, [QSI(0,Fr(-3,8),0,Fr(-1,8)), QSI(0,0,0,Fr(-3,8)), QSI(0,Fr(1,4),0,0)],
         [QSI(0,Fr(-3,8),0,Fr(1,8)), QSI(0,0,0,Fr(3,8)), QSI(0,Fr(1,4),0,0)], QSI(Fr(1,8),Fr(-1,8))]
s30 = [{}, [QSI(0,Fr(1,12))], [QSI(0,Fr(1,12))], QSI(Fr(-1,12))]
s21 = [{}, [QSI(0,Fr(-1,4))], [QSI(0,Fr(-1,4))], QSI(Fr(1,4))]

# sanity: all vanish at v=b (Gm property)
for nm, G in [('G3',G3),('G21',G21),('s30',s30),('s21',s21)]:
    print(nm, "Gm(b) =", pf_eval_at(G, BQ))

def simple_of(P):
    s, h, c = pf_split(P)
    assert c.is_zero()
    return s, h

def show(P, label):
    real = {str(r.sympy()): str(c.sympy()) for r, c in P[0].items()}
    ip = [str(e.sympy()) for e in P[1]]
    print(f"  {label}: real={real}  ipos(order1..)={ip}")

# C1: SimplePart[-3 G3 DL1 - G21 DL2] == R_20(3)
L = pf_add(pf_scale(pf_mul_DL(G3, DL1), Fr(-3)), pf_neg(pf_mul_DL(G21, DL2)))
s20, h20 = simple_of(L)
print("C1 LHS simple:"); show(s20, "simple(-3G3 DL1 - G21 DL2)")
R20_3 = pf_from_poly(*__import__('buildintegrand').__dict__.get('x', ([],[]))) if False else None
# build expected R_20(3) as PF
def pf_real(pairs):
    P = pf_zero()
    for r, c in pairs:
        P[0][r] = c
    return P
R20_3exp = pf_real([(ZERO, QSI(Fr(1,4),Fr(-1,4))), (AQ, QSI(Fr(1,4))), (-AQ, QSI(Fr(-1,4))), (-BQ, QSI(Fr(-1,2)))])
print("C1:", "PASS" if not (pf_sub(s20, R20_3exp)[0] or pf_sub(s20,R20_3exp)[1] or pf_sub(s20,R20_3exp)[2]) else "FAIL")

# C2: SimplePart[-3 s30 DL1 - s21 DL2] == S_20 = 1/(2(v+a)) - 1/(2(v+b))
L2 = pf_add(pf_scale(pf_mul_DL(s30, DL1), Fr(-3)), pf_neg(pf_mul_DL(s21, DL2)))
sS20, hS20 = simple_of(L2)
S20exp = pf_real([(-AQ, QSI(Fr(1,2))), (-BQ, QSI(Fr(-1,2)))])
d = pf_sub(sS20, S20exp)
print("C2:", "PASS" if not (d[0] or d[1] or d[2]) else "FAIL"); show(sS20, "simple slope20")

# C3: SimplePart[-2 G21 DL1] == R_11(3) = (1+sqrt2)/(4v) - 1/(2(v+a))
L3 = pf_scale(pf_mul_DL(G21, DL1), Fr(-2))
s11, h11 = simple_of(L3)
R11_3exp = pf_real([(ZERO, QSI(Fr(1,4),Fr(1,4))), (-AQ, QSI(Fr(-1,2)))])
d = pf_sub(s11, R11_3exp)
print("C3:", "PASS" if not (d[0] or d[1] or d[2]) else "FAIL"); show(s11, "simple(-2G21 DL1)")

# C4: SimplePart[-2 s21 DL1] == S_11 = 1/(2v) - 1/(v+a)
L4 = pf_scale(pf_mul_DL(s21, DL1), Fr(-2))
sS11, hS11 = simple_of(L4)
S11exp = pf_real([(ZERO, QSI(Fr(1,2))), (-AQ, QSI(-1))])
d = pf_sub(sS11, S11exp)
print("C4:", "PASS" if not (d[0] or d[1] or d[2]) else "FAIL"); show(sS11, "simple slope11")
