# 07_area_checks.py -- validate the exact 1-D area formula (superlib.py) against
#  (i) exact values p=1 (octahedron, S=4 sqrt3) and p=2 (sphere, S=4 pi),
#  (ii) a completely independent direct surface integral (graph z=(1-x^p-y^p)^(1/p) over the
#       region where z is the largest coordinate, times 3 and symmetrised), p = 3, 4, 10,
#  (iii) the squircle perimeter P(p)=8-4 eps F(eps) against direct arc length,
#  (iv) the series for h_eps against direct quadrature.
import mpmath as mp, sys
from superlib import h, H, A_oct, S_area, F2, Qfun
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
D = mp.mp.dps

print("(iv) h_eps(b): series/asymptotic vs quadrature")
for eps in ['0.3', '0.05', '-0.2']:
    for b in ['0.2', '1.5', '6']:
        e = mp.mpf(eps); bb = mp.mpf(b)
        hq = mp.quad(lambda x: -mp.expm1(-x**(2-2*e)/bb**2)*x**(e-1)*mp.exp(-x), [0, bb**(1/(1-e)), mp.inf])
        print("  eps=%5s b=%4s  h=%s  diff=%s" % (eps, b, mp.nstr(h(e, bb), 20), mp.nstr(h(e, bb)-hq, 3)))

print("(i) exact cases")
S1 = S_area(1); S2 = S_area(2)
print("  p=1: S =", mp.nstr(S1, D), " 4sqrt3 =", mp.nstr(4*mp.sqrt(3), D), " diff", mp.nstr(S1-4*mp.sqrt(3), 3))
print("  p=2: S =", mp.nstr(S2, D), " 4pi    =", mp.nstr(4*mp.pi, D), " diff", mp.nstr(S2-4*mp.pi, 3))

def S_direct(p):
    p = mp.mpf(p)
    def f(x, y):
        zp = 1 - x**p - y**p
        rx = (x**p/zp)**(2-2/p); ry = (y**p/zp)**(2-2/p)
        return mp.sqrt(1 + rx + ry)
    x3 = 3**(-1/p); x2 = 2**(-1/p)
    def inner(x):
        ym = x if x <= x3 else (1-2*x**p)**(1/p)
        pts = [0] + [ym*(1-mp.mpf(1)/2**k) for k in range(1, 8)] + [ym]
        return mp.quad(lambda y: f(x, y), pts)
    xpts = [0] + [x3*(1-mp.mpf(1)/2**k) for k in range(1, 8)] + [x3] + [x3+(x2-x3)*(1-mp.mpf(1)/2**k) for k in range(1, 6)] + [x2]
    return 48*mp.quad(inner, xpts)

print("(ii) independent direct surface integral")
with mp.workdps(20):
    for p in [3, 4, 10]:
        Sd = mp.re(S_direct(p))
        Sf = S_area(p)
        print("  p=%d  direct=%s  formula=%s  diff=%s" % (p, mp.nstr(Sd, 18), mp.nstr(Sf, 18), mp.nstr(Sd-Sf, 3)))

print("(iii) squircle perimeter check P = 8 - 4 eps F(eps)")
for p in [3, 10, 50]:
    pp = mp.mpf(p); e = 1/pp
    # independent: 8 x arc length of the graph y=(1-x^p)^(1/p) for 0<=x<=2^(-1/p) (slope <= 1)
    yp = lambda x: -x**(pp-1)*(1-x**pp)**(1/pp-1)
    x2 = 2**(-1/pp)
    Pd = 8*mp.quad(lambda x: mp.sqrt(1+yp(x)**2), [0] + [x2*(1-mp.mpf(1)/2**k) for k in range(1, 10)] + [x2])
    Pf = 8 - 4*e*F2(e)
    print("  p=%d  direct=%s  formula=%s  diff=%s" % (p, mp.nstr(Pd, D), mp.nstr(Pf, D), mp.nstr(Pd-Pf, 3)))
