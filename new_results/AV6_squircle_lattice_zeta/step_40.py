code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, quad, exp, log, power, gamma
import time
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
ap = gamma(1+1/p); ga = gamma(a)
cache = {}
def mpw(m):
    v = cache.get(m)
    if v is None:
        v = power(m, p); cache[m] = v
    return v
def theta(t):
    f = mpf(0); cutoff = mpf(10)**(-85); m = 0
    while True:
        m += 1
        e = exp(-t*mpw(m)); f += e
        if e < cutoff: break
    return 1+2*f
def F(t):
    th = theta(t); return th*th-1
def total(delta):
    ld = float(log(delta))
    small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
    pts = [ld]
    while pts[-1] < -25: pts.append(pts[-1]+25)
    pts.append(0.0)
    I1 = mpf(0)
    for i in range(len(pts)-1):
        I1 += quad(lambda x: exp(a*x)*F(exp(x)), [pts[i], pts[i+1]], method='gauss-legendre')
    I2 = quad(lambda x: exp(a*x)*F(exp(x)), [0.0, 5.0])
    return small + I1 + I2
t0=time.time()
ZA_ = ZA(50, 4)
print("ZA*ga   :", nstr(ZA_*ga, 45))
kap = mpf('0.189179')
for dexp in [150, 120, 100]:
    dd = mpf(10)**(-dexp)
    rex = kap*dd**(-1/mpf(49))
    T = total(dd)
    print(f"delta=1e-{dexp} (R~e^-{nstr(rex,3)}): {nstr(T, 45)}   ({time.time()-t0:.0f}s)")
'''
exec(code)
