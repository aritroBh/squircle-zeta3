import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
def r_oct(phi, p):
    c, s = cos(phi), sin(phi)
    return (c**p + s**p + mpf(2)**(-p/2)*((c+s)**p + (c-s)**p))**(-1/p)
def red(alpha):
    # distance to nearest multiple of pi/4
    x = (alpha + pi/8) % (pi/4)
    return abs(x - pi/8)
def h_oct(theta, p):
    def hof(phi):
        return r_oct(phi, p)*cos(red(theta-phi))
    lo, hi = mpf(0), pi/4
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = hof(x1), hof(x2)
    for _ in range(120):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = hof(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = hof(x1)
    return (f1+f2)/2
print("p=2:", nstr(h_oct(mpf('0.6'),mpf(2)),25), nstr(1/sqrt(2),25))
print("p=2:", nstr(h_oct(mpf('1.1'),mpf(2)),25))
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(th,p), [0, pi/16, pi/8])
t0=time.time()
print("P(2):", nstr(P_oct(mpf(2)),30), nstr(pi*sqrt(2),30))
print("time", time.time()-t0)
