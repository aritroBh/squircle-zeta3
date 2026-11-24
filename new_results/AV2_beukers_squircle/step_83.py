import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
from mpmath import mp, mpf, quad, log, log1p, exp, nstr
mp.dps = 100
exec(open('av2core.py').read())
cs2 = shifted_legendre_coeffs(2)
cs3 = shifted_legendre_coeffs(3)
def P(x, cs):
    return sum(c*x**i for i, c in enumerate(cs))
def L_direct(n, p):
    cs = shifted_legendre_coeffs(n)
    p = mpf(p); ep = 1/p
    def f(s, t):
        lA = ep*log1p(-s**p)
        A = exp(lA)
        y = t*A
        xy = s*y
        return A * (-log(xy))/(1-xy) * P(s, cs)*P(y, cs)
    return quad(f, [0, 1], [0, 1])
Mdict = {}
for j in range(13):
    for k in range(j, 13):
        Mdict[(j,k)] = M_moment(j, k, 2, mmin=500)
print("L_2(2) moments:", nstr(L_legendre(2, 2, Mdict), 25))
print("L_2(2) direct :", nstr(L_direct(2, 2), 25))
print("L_3(2) moments:", nstr(L_legendre(3, 2, Mdict), 25))
print("L_3(2) direct :", nstr(L_direct(3, 2), 25))
