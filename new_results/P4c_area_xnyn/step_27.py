# step_25.py : rho residual trend + symmetry check E_ρ(n) = E_{1/ρ}(ρn)
import math, numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi, zeta
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
    sgmax = log(n) + 50
    npan = int(math.ceil(float(sgmax)/4.0))
    total = mpf(0)
    for j in range(npan):
        s0 = 4*mpf(j); s1 = min(4*mpf(j+1), sgmax)
        if s0 >= sgmax: break
        def outer(sg):
            def gfun(d):
                u = (sg+d)/2; w = (sg-d)/2
                W = n*exp(-sg)*sqrt(exp(2*e*u) + rho**2*exp(2*e*w/rho))
                return (sqrt(1+W*W)-1)*exp(-e*u)*exp(-e*w/rho)*mpf('0.5')
            return gli(gfun, -sg, sg, deg=20)
        total += gli(outer, s0, s1, deg=20)
    return total/(rho*n**2)
# symmetry check
E2_1000 = E_rho_direct(1000, 2)
Ehalf_2000 = E_rho_direct(2000, mpf('0.5'))
print("E_2(1000)       =", nstr(E2_1000, 30))
print("E_{1/2}(2000)   =", nstr(Ehalf_2000, 30), "  diff =", nstr(E2_1000-Ehalf_2000, 8))
# trend of R_ρ diffs for rho=0.5 and rho=0.25
L2 = log(2); z2 = pi**2/6
def R_rho_pred(n, rho):
    rho = mpf(rho); L = log(mpf(n))
    lam = log(1+rho**2)/2
    c1 = -(1+rho**3)/(rho*sqrt(1+rho**2))
    mt1 = 1-L2; mt2 = 1-L2+L2**2/2+z2/2
    return c1/rho - ((L+lam)**2/2 + mt1*(L+lam) + mt2)/rho
for rho in ['0.5','0.25']:
    print("rho =", rho)
    for n in [2000, 4000, 8000, 16000]:
        Ed = E_rho_direct(n, rho)
        Rnum = n*(n*Ed - sqrt(1+mpf(rho)**2)/mpf(rho))
        diff = Rnum - R_rho_pred(n, rho)
        print("  n=%5d: diff=%s   (log^3 n)/n=%s   diff/(that)=%s" % (n, nstr(diff,12), nstr(log(n)**3/n,8), nstr(diff/(log(n)**3/n),8)))
