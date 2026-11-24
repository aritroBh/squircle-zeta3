# step_65: stability of L3 at extreme rational arguments
import mpmath as mp
from fractions import Fraction as Fr
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def test(z):
    for dps in (60, 100):
        mp.mp.dps = dps
        def Li3(x): return mp.re(mp.polylog(3, x))
        def Li2(x): return mp.re(mp.polylog(2, x))
        def L3mpf(q):
            if q > 1: return L3mpf(1/q)
            return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
        print(f"  L3({z}) dps={dps}: {mp.nstr(L3mpf(mpfq(z)), 30)}")
for z in ['-4374','-2400','-224','-125','36','28','21','2']:
    test(z)
