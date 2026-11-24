# step_21.py : z = x^n y^{ρn} : verify E = nJ_ρ - C_ρ and the R_ρ asymptotic formula
import math, numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi, zeta, binomial
mp.dps = 60
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

def E_rho_direct(n, rho):
    n = mpf(n); rho = mpf(rho); e = 1/n
    sgmax = log(n) + 60
    npan = int(math.ceil(float(sgmax)/4.0))
    total = mpf(0)
    for j in range(npan):
        s0 = 4*mpf(j); s1 = min(4*mpf(j+1), sgmax)
        if s0 >= sgmax: break
        def outer(sg):
            def g(d):
                u = (sg+d)/2; w = (sg-d)/2
                W = n*exp(-sg)*sqrt(exp(2*e*u) + rho**2*exp(2*e*w/rho))
                return (sqrt(1+W*W)-1)*exp(-e*u)*exp(-e*w/rho)*mpf('0.5')
            return gli(g, -sg, sg, deg=24)
        total += gli(outer, s0, s1, deg=24)
    return total/(rho*n**2)

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)
def C_rho(n, rho):
    n = mpf(n); rho = mpf(rho); e = 1/n
    sgmax = log(n) + 60
    npan = int(math.ceil(float(sgmax)/4.0))
    total = mpf(0)
    for j in range(npan):
        s0 = 4*mpf(j); s1 = min(4*mpf(j+1), sgmax)
        if s0 >= sgmax: break
        def outer(sg):
            def g(d):
                u = (sg+d)/2; w = (sg-d)/2
                W = n*exp(-sg)*sqrt(exp(2*e*u) + rho**2*exp(2*e*w/rho))
                return fW(W)*exp(-e*u)*exp(-e*w/rho)*mpf('0.5')
            return gli(g, -sg, sg, deg=20)
        total += gli(outer, s0, s1, deg=24)
    return total/(rho*n**2)

def J_series(n, rho, kmax=300):
    n = mpf(n); rho = mpf(rho); m = rho*n
    s = mpf(0)
    for k in range(kmax):
        b = binomial(mpf('0.5'), k)
        t1 = rho**(2*k)/(m-2*k+1)*(1/(n+2*k) - rho**(m-2*k+1)/(m+n+1))
        t2 = rho**(m+1)/((m+2*k)*(m+n+1))
        s += b*(t1+t2)
    return s

L2 = log(2); z2 = pi**2/6
mt1 = 1-L2; mt2 = 1-L2 + L2**2/2 + z2/2
def R_rho(n, rho):
    rho = mpf(rho); L = log(mpf(n))
    lam = log(1+rho**2)/2
    c1 = -(1+rho**3)/(rho*sqrt(1+rho**2))
    return c1/rho - (L**2/2 + mt1*L + mt2 + lam*L + lam*mt1 + lam**2/2)/rho
    # = c1/ρ − [½(L+λ)² + m̃1(L+λ) + m̃2]/ρ

for rho in ['0.5','2']:
    for n in [1000, 10000]:
        Ed = E_rho_direct(n, rho)
        Eid = n*J_series(n, rho) - C_rho(n, rho)
        lead = sqrt(1+mpf(rho)**2)/(mpf(rho)*n)
        Rnum = n*(n*Ed - sqrt(1+mpf(rho)**2)/mpf(rho))
        Rpred = R_rho(n, rho)
        print("rho=%s n=%d:" % (rho, n))
        print("  E direct      =", nstr(Ed, 30))
        print("  nJ - C_rho    =", nstr(Eid, 30), "  diff =", nstr(Ed-Eid, 8))
        print("  leading term  =", nstr(lead, 30))
        print("  R_rho num     =", nstr(Rnum, 25))
        print("  R_rho pred    =", nstr(Rpred, 25), "  diff =", nstr(Rnum-Rpred, 10), "  ~O(log^3 n)/n =", nstr(log(n)**3/n, 6))
