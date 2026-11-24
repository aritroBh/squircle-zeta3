
import json
from mpmath import mp, mpf, zeta, fabs, nstr, log
mp.dps = 100
exec(open('av2core.py').read())
from math import comb, gcd
from functools import reduce
def lcm(a, b): return a*b//gcd(a, b)
def L(n): return reduce(lcm, range(1, n+1), 1)
def apery(n): return sum(comb(n,k)**2 * comb(n+k,k)**2 for k in range(n+1))

# cross-check L_2(2) by direct 2D quadrature
print("=== direct 2D check of L_2(2) ===")
cs = shifted_legendre_coeffs(2)
def P(x, cs):
    return sum(c*x**i for i, c in enumerate(cs))
def L2_direct(p):
    p = mpf(p); ep = 1/p
    def f(s, t):
        lA = ep*log1p(-s**p); lB = ep*log1p(-t**p)
        A = exp(lA); B = exp(lB)
        xy = s*t*A
        # -log(xy)/(1-xy): xy > 0 interior
        return A * (-log(xy))/(1-xy) * P(s, cs)*P(t*B, cs)
    from mpmath import quad
    return quad(f, [0, 1], [0, 1])
Mdict = {}
for j in range(13):
    for k in range(j, 13):
        Mdict[(j,k)] = M_moment(j, k, 2, mmin=500)
print("L_2(2) via moments:", nstr(L_legendre(2, 2, Mdict), 30))
print("L_2(2) direct      :", nstr(L2_direct(2), 30))

results = {}
for p in [1, 2, 4, 8, 32]:
    mm = int(220*p)+50
    Mdict = {}
    for j in range(13):
        for k in range(j, 13):
            Mdict[(j,k)] = M_moment(j, k, p, mmin=mm)
    B2p, _ = B2_series(p)
    print("=== p = %s (mmin=%d) ===" % (p, mm))
    rows = []
    for n in range(0, 13):
        Ln = L_legendre(n, p, Mdict)
        an = apery(n)
        r_n = Ln - 2*an*B2p
        rate = fabs(Ln)**(1.0/n) if n else mpf(0)
        print("n=%2d L_n=%s rate=%s r_n=%s" % (n, nstr(Ln, 18), nstr(rate, 9), nstr(r_n, 18)))
        rows.append(dict(n=n, Ln=mp.nstr(Ln, 55), rn=mp.nstr(r_n, 55), rate=mp.nstr(rate, 16)))
    results[str(p)] = dict(B2=mp.nstr(B2p, 60), rows=rows)
    json.dump(results, open('legendre_finite.json','w'), indent=1)
print("saved")
