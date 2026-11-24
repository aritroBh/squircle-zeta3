
import os, json, sys
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
sys.path.insert(0, os.path.join(BASE,'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
from qsqrt2 import QS2
from polylog_exact import P, LAM, LL, PI2, Z3, LI3, LI2
import mpmath as mp
mp.mp.dps = 70
B = QS2(-1,1)

# rebuild W3(2) from the printed form
W3_2 = (Z3 * QS2(Fr(-19,16), Fr(189,64)) + LI3(-B) * QS2(11,-15) + LI3(B) * QS2(11,15)
        + LI3(QS2(0,Fr(1,2))) * QS2(0,-27) + PI2*LAM * QS2(Fr(-1,24), Fr(-65,48))
        + LAM*LI2(B) * QS2(0, Fr(-7,2)) + LL*PI2 * QS2(Fr(-7,12), Fr(15,8)) + LL*LI2(B) * 12
        + LAM*LAM*LAM * QS2(Fr(1,6), Fr(9,16)) + LAM*LAM*LL * QS2(0, Fr(11,8))
        + LAM*LL*LL * QS2(-1, Fr(-7,4)) + LL*LL*LL * QS2(Fr(25,6), Fr(-7,3)))
W3_M2 = (Z3 * QS2(Fr(-19,16), Fr(189,64)) + LI3(B) * QS2(11, 15) + LI3(-B) * QS2(11, -15)
         + LI3(QS2(0, Fr(1,2))) * QS2(0, -27) + LI2(B) * (LL * 12 - LAM * QS2(0, Fr(7,2)))
         + PI2 * LL * QS2(Fr(-7,12), Fr(15,8)) + PI2 * LAM * QS2(Fr(-1,24), Fr(-65,48))
         + LAM*LAM*LAM * QS2(Fr(1,6), Fr(9,16)) + LAM*LAM*LL * QS2(0, Fr(11,8))
         + LAM*LL*LL * QS2(-1, Fr(-7,4)) + LL*LL*LL * QS2(Fr(25,6), Fr(-7,3)))
print("W3(2) == thm:M2 weight-3 part:", "PASS" if not (W3_2 - W3_M2).d else "FAIL")

vals = json.load(open('w3_values.json'))['values']
def p_from_json(lst):
    out = P()
    for atoms, c in lst:
        term = P.const(QS2(Fr(c[0]), Fr(c[1])))
        for a in atoms:
            if isinstance(a, str):
                term = term * P({(a,): 1})
            else:
                z = QS2(Fr(a[1][0]), Fr(a[1][1]))
                term = term * (LI2(z) if a[0]=='Li2' else LI3(z))
        out = out + term
    return out
W3_3 = p_from_json(vals['3'])
DELTA2 = W3_3 - W3_2
print("W3(3) - W3(2) =", DELTA2)
print("numeric:", mp.nstr(DELTA2.num(), 30))
print("cf k>=3 step DELTA = -1/8*z3 - 1/24*lam*pi2 - 1/12*lam^3 =", mp.nstr((Z3*Fr(-1,8) - PI2*LAM*Fr(1,24) - LAM*LAM*LAM*Fr(1,12)).num(), 30))
