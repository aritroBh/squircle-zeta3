import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr
mp.dps = 80
def h_lame(theta, p, lam):
    if theta == 0 or theta == pi/2: return mpf(1)
    c, s = cos(theta), sin(theta)
    def F(x):
        return c - (1/lam)*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*s
    lo, hi = mpf(0), 1-mpf(10)**(-60)
    for _ in range(260):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*c + (1-x**p)**(1/(lam*p))*s
p, lam = mpf(20), mpf(2)
for splits in [[0, pi/2], [0, pi/8, pi/4, 3*pi/8, pi/2], [0, pi/16, pi/8, 3*pi/16, pi/4, 5*pi/16, 3*pi/8, 7*pi/16, pi/2]]:
    P = 4*mpmath.quad(lambda th: h_lame(th, p, lam), splits)
    print(nstr(P, 60))
