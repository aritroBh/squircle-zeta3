"""step_32.py -- filter relations with row[0] != 0; verify ladder by differentiation."""
import mpmath as mp
from ir_utils import find_relations
from fractions import Fraction
mp.mp.dps = 100

def Hval(z, xi):
    return mp.quad(lambda w: mp.polylog(2, w)/(w - xi), [0, z/2, z])

def candidates(z, xi):
    L = mp.log
    return [
     ("Li3(z/xi)", mp.polylog(3, z/xi)), ("Li3((xi-z)/xi)", mp.polylog(3, (xi-z)/xi)),
     ("Li3(z/(z-xi))", mp.polylog(3, z/(z-xi))), ("Li3(xi/(xi-z))", mp.polylog(3, xi/(xi-z))),
     ("LzL1zxi2", L(z)*L(1-z/xi)**2), ("L1zxi3", L(1-z/xi)**3),
     ("L1zxi2Lxi", L(1-z/xi)**2*L(xi)), ("zeta3", mp.zeta(3)),
    ]

xi = mp.mpf('2.7'); z = mp.mpf('0.31')
Hz = Hval(z, xi)
cand = candidates(z, xi)
vals = [mp.re(v) for _, v in cand]
rels = find_relations([Hz] + vals, maxcoeff=10**5, resid_decades=20, max_report=4)
formula = None
for row, mc, res in rels:
    if row[0] == 0:
        print("internal candidate relation (not involving H)")
        continue
    q = [Fraction(-row[i], row[0]) for i in range(1, len(row))]
    formula = q
    print("H-relation resid", mp.nstr(res, 2), " H =", [(cand[i][0], str(q[i])) for i in range(len(q)) if q[i] != 0])
    break

def Hcand(z, xi):
    cand = candidates(z, xi)
    return mp.fsum(mp.mpf(q.numerator)/q.denominator*mp.re(v) for q, (_, v) in zip(formula, cand))

print("\ndifferentiation check:")
import random
random.seed(2)
worst = 0
for _ in range(8):
    zz = mp.mpf(random.uniform(0.05, 0.95))
    xx = mp.mpf(random.uniform(1.2, 4.0))
    h = mp.mpf(10)**(-40)
    der = (Hcand(zz + h, xx) - Hcand(zz - h, xx))/(2*h)
    tgt = mp.polylog(2, zz)/(zz - xx)
    worst = max(worst, abs(der - tgt))
    print(f"  z={mp.nstr(zz,6)}, xi={mp.nstr(xx,6)}: |der-target| = {mp.nstr(abs(der-tgt), 2)}")
print("value check:", mp.nstr(Hcand(mp.mpf('0.31'), mp.mpf('2.7')) - Hz, 3))
