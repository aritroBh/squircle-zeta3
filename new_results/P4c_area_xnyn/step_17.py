# step_15.py : high-precision E(n) for n=1e3..1e7 via exact decomposition (rerun of prev step_10)
import time, math
import numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi
mp.dps = 70
_gl_cache = {}
def gl_nodes(deg):
    if deg in _gl_cache: return _gl_cache[deg]
    xs, ws = np.polynomial.legendre.leggauss(deg)
    nodes = []
    for x0 in xs:
        x = mpf(float(x0))
        for _ in range(60):
            p0, p1 = mpf(1), x
            for k in range(2, deg+1):
                p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
            dP = deg*(p0 - x*p1)/(1-x*x)
            dx = p1/dP; x -= dx
            if abs(dx) < mpf(10)**(-mp.dps+10): break
        p0, p1 = mpf(1), x
        for k in range(2, deg+1):
            p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
        dP = deg*(p0 - x*p1)/(1-x*x)
        nodes.append((x, 2/((1-x*x)*dP*dP)))
    _gl_cache[deg] = nodes
    return nodes
def gli(func, a, b, deg):
    mid=(a+b)/2; half=(b-a)/2
    return sum(w*func(mid+half*x) for x,w in gl_nodes(deg))*half
def Phi(n):
    n = mpf(n); width = mpf(10)/n; total = mpf(0)
    for j in range(16):
        total += gli(lambda y: exp(-n*y)*sqrt(1+exp(-2*y)), j*width, (j+1)*width, deg=24)
    return total
def C_int(n, taildigits=60):
    n = mpf(n); e = 1/n; aa = 1 - e/2
    def f_w(w): return 1 + w - sqrt(1+w*w)
    def inner(sg):
        base = n*exp(-aa*sg)
        return gli(lambda d: f_w(base*sqrt(2*cosh(e*d))), 0, sg, deg=22)
    L = log(n); ln10 = log(10)
    sgmax = L + (taildigits+6)*ln10
    npan = int(math.ceil(float(sgmax)/3.0))
    total = mpf(0)
    for j in range(npan):
        a0 = 3*mpf(j); b0 = min(3*mpf(j+1), sgmax)
        if a0 >= sgmax: break
        total += gli(lambda sg: exp(-e*sg)*inner(sg), a0, b0, deg=30)
    return total

L2 = log(2); z2 = zeta(2) if False else pi**2/6
beta = L2/2 - 1
gamma_ = -(1+sqrt(2)) + L2/2 - pi**2/12 - L2**2/8
def R(n):
    ln = log(mpf(n))
    return -ln**2/2 + beta*ln + gamma_
results = {}
t0 = time.time()
for k in [3,4,5,6,7]:
    n = 10**k
    C = C_int(n)
    E = Phi(n)/(1+1/(2*mpf(n))) - C/mpf(n)**2
    resid = mpf(n)**2*(E - sqrt(2)/n) - R(n)
    results[n] = (E, resid, C)
    print("n=1e%d E = %s" % (k, nstr(E, 55)))
    print("   resid n^2(E-sqrt2/n)-R(n) = %s  (%.1fs)" % (nstr(resid, 30), time.time()-t0))
import pickle
pickle.dump({k:(str(v[0]),str(v[1]),str(v[2])) for k,v in results.items()}, open("E_values_new.pkl","wb"))
