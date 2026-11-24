from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler)
mp.dps = 100
def N(x,d=30): return nstr(x,d)
Lp = log(pi)
eta4 = mpf(7)/8*zeta(4)
etap4 = (log(2)/8)*zeta(4) + (mpf(7)/8)*zeta(4, derivative=1)
A3_direct = nsum(lambda n: (-1)**n*(ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n**4, [1, inf])
CiS4 = nsum(lambda n: (-1)**n*ci(n*pi)/n**4, [1, inf])
A3_pieces = CiS4 + euler*eta4 - Lp*zeta(4) + etap4
print("A3 direct =", N(A3_direct,38))
print("A3 pieces =", N(A3_pieces,38), " diff:", nstr(abs(A3_direct-A3_pieces),3))
# sub-checks
print("sum (-1)^n log n /n^4 =", N(nsum(lambda n: (-1)**n*log(n)/n**4, [1,inf]),30), " vs -eta'(4) =", N(-etap4,30))
print("sum (-1)^n ci /n^4 partial:", N(nsum(lambda n: (-1)**n*ci(n*pi)/n**4, [1,20000]),20))
