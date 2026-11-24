import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr
mp.dps = 80
def h_lame(theta, p, lam):
    if theta == 0 or theta == pi/2: return mpf(1)
    c, s = cos(theta), sin(theta)
    inv_lam = 1/lam
    def F(x):
        return c - inv_lam*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*s
    lo, hi = mpf(0), 1-mpf(10)**(-60)
    for _ in range(260):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*c + (1-x**p)**(1/(lam*p))*s
def P_lame_sup(p, lam):
    return 4*mpmath.quad(lambda th: h_lame(th, p, lam), [0, pi/8, pi/4, 3*pi/8, pi/2])
def P_lame_arc(p, lam):
    a = 2/p; b = 2/(lam*p)
    def f(t):
        dx = -a*cos(t)**(a-1)*sin(t)
        dy = b*sin(t)**(b-1)*cos(t)
        return sqrt(dx**2+dy**2)
    return 4*mpmath.quad(f, [0, pi/4, pi/2])
for p,lam in [(mpf(5),mpf(2)), (mpf(3),mpf('0.5')), (mpf(20),mpf(2)), (mpf(100),mpf(2))]:
    Pa = P_lame_sup(p,lam); Pb = P_lame_arc(p,lam)
    print(f"p={p}, lam={lam}: support {nstr(Pa,45)}")
    print(f"            arclength {nstr(Pb,45)}  diff {nstr(Pa-Pb,4)}")
