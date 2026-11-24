"""Planner check of ZQ5's D=5 formula, own code: Li5 by direct series, L(chi5,5) by Hurwitz zeta at two precisions.
L(chi5,5) = 5^-5 [zeta(5,1/5) - zeta(5,2/5) - zeta(5,3/5) + zeta(5,4/5)]  (chi5 = Legendre symbol mod 5).
Exact symbol check: one letter (log phi) plus log 2; done by hand-coded integer arithmetic below."""
import sys
from fractions import Fraction as F
import mpmath as mp
def li5(x, eps):
    t, xn, n = mp.mpf(0), mp.mpf(1), 1
    while True:
        xn *= x; t += xn / mp.mpf(n) ** 5
        if abs(xn) < eps: return t
        n += 1
def run(d):
    mp.mp.dps = d + 30; eps = mp.mpf(10) ** (-(d + 25))
    phi = (1 + mp.sqrt(5)) / 2; lp = mp.log(phi); r = 1 / phi
    rhs = mp.sqrt(5) * (F(1728, 6875) * 0 + mp.mpf(1728) / 6875 * (li5(r, eps) - li5(-r, eps))
          - mp.mpf(512) / 6875 * li5(r ** 3, eps) + mp.mpf(16) / 6875 * li5(r ** 6, eps)
          + mp.mpf(648) / 34375 * lp ** 5 - mp.mpf(48) / 6875 * mp.pi ** 2 * lp ** 3 + mp.mpf(64) / 20625 * mp.pi ** 4 * lp)
    L = (mp.zeta(5, mp.mpf(1) / 5) - mp.zeta(5, mp.mpf(2) / 5) - mp.zeta(5, mp.mpf(3) / 5) + mp.zeta(5, mp.mpf(4) / 5)) / mp.mpf(5) ** 5
    return L, rhs
# exact symbol (coefficients of phi^{x5} and of 2 (x) phi^{x4}); v(x), v(1-x) in (log phi, log 2) from x = phi^-k
terms = [(F(1728, 6875), (-1, 0), (-2, 0)), (F(-1728, 6875), (-1, 0), (1, 0)),
         (F(-512, 6875), (-3, 0), (-2, 1)), (F(16, 6875), (-6, 0), (-3, 2))]
s_phi5 = sum(c * (-v1[0]) * v[0] ** 4 for c, v, v1 in terms) + F(648, 34375) * 120
s_2phi4 = sum(c * (-v1[1]) * v[0] ** 4 for c, v, v1 in terms)
print('exact symbol: phi^5 coeff', s_phi5, ' 2(x)phi^4 coeff', s_2phi4)
for d in (int(sys.argv[1]) if len(sys.argv) > 1 else 1500, 600):
    L, rhs = run(d)
    print(f'dps {d}: L(chi5,5) = {mp.nstr(L, 30)}  RHS - L = {mp.nstr(rhs - L, 6)}')
