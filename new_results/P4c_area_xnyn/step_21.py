# step_19.py : is C_int(1e6) accurate? recompute with tighter params; also C0 sigma-form at 1e6
import math, numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi
mp.dps = 80
_gl_cache = {}
def gl_nodes(deg):
    if deg in _gl_cache: return _gl_cache[deg]
    xs, ws = np.polynomial.legendre.leggauss(deg)
    nodes = []
    for x0 in xs:
        x = mpf(float(x0))
        for _ in range(80):
            p0, p1 = mpf(1), x
            for k in range(2, deg+1):
                p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
            dP = deg*(p0 - x*p1)/(1-x*x); dx = p1/dP; x -= dx
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
def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)

def C_int(n, taildigits=65, inner_deg=36, outer_deg=44, panelw=2):
    n = mpf(n); e = 1/n; aa = 1 - e/2
    def inner(sg):
        base = n*exp(-aa*sg)
        return gli(lambda d: fW(base*sqrt(2*cosh(e*d))), 0, sg, deg=inner_deg)
    L = log(n); ln10 = log(10)
    sgmax = L + (taildigits+6)*ln10
    npan = int(math.ceil(float(sgmax)/panelw))
    total = mpf(0)
    for j in range(npan):
        a0 = panelw*mpf(j); b0 = min(panelw*mpf(j+1), sgmax)
        if a0 >= sgmax: break
        total += gli(lambda sg: exp(-e*sg)*inner(sg), a0, b0, deg=outer_deg)
    return total

def C0_sig(n):
    n = mpf(n)
    sgmax = log(n*sqrt(2)) + 50
    npan = int(math.ceil(float(sgmax)/3.0))
    return sum(gli(lambda sg: sg*fW(n*sqrt(2)*exp(-sg)), 3*j, min(3*(j+1), sgmax), 40) for j in range(npan) if 3*j < sgmax)

import time
n = mpf(10)**6
t0 = time.time()
C0 = C0_sig(n)
print("C0(1e6) sigma-form =", nstr(C0, 50), "(%.1fs)"%(time.time()-t0))
L2 = log(2); z2 = pi**2/6
mt1 = 1-L2; mt2 = 1-L2 + L2**2/2+z2/2
L = log(n*sqrt(2))
poly1 = L**2/2 + mt1*L + mt2
print("poly1 - 1/(2N)     =", nstr(poly1 - 1/(2*n*sqrt(2)), 50))
print("C0 - (poly1-1/2N)  =", nstr(C0 - (poly1 - 1/(2*n*sqrt(2))), 15))
t0 = time.time()
Ctight = C_int(n)
print("C_int(1e6) tight   =", nstr(Ctight, 50), "(%.1fs)"%(time.time()-t0))
print("prev C_int(1e6)    = 105.9965931431501241549511127430218278965")
print("diff               =", nstr(Ctight - mpf('105.9965931431501241549511127430218278965'), 15))
