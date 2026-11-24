# step_24.py : final validation batch: planner cross-check, Phi Watson, A/B constants, rho at 0.25 and 3
import math, numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi, zeta, binomial, inf, quad
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
def Phi(n):
    n = mpf(n); width = mpf(10)/n; total = mpf(0)
    for j in range(16):
        total += gli(lambda y: exp(-n*y)*sqrt(1+exp(-2*y)), j*width, (j+1)*width, deg=24)
    return total
def C_int(n, taildigits=55):
    n = mpf(n); e = 1/n; aa = 1 - e/2
    def f_w(w): return 1 + w - sqrt(1+w*w)
    def inner(sg):
        base = n*exp(-aa*sg)
        return gli(lambda d: f_w(base*sqrt(2*cosh(e*d))), 0, sg, deg=22)
    L = log(n); sgmax = L + (taildigits+6)*log(10)
    npan = int(math.ceil(float(sgmax)/3.0))
    total = mpf(0)
    for j in range(npan):
        a0 = 3*mpf(j); b0 = min(3*mpf(j+1), sgmax)
        if a0 >= sgmax: break
        total += gli(lambda sg: exp(-e*sg)*inner(sg), a0, b0, deg=30)
    return total
# 1) planner cross-check E(10), E(100)
for n in [10, 100]:
    E = Phi(n)/(1+1/(2*mpf(n))) - C_int(n)/mpf(n)**2
    print("n=%d: nE = %s" % (n, nstr(n*E, 40)))
print("planner: nE(10) = 0.81901914523902943304880 ; nE(100) = 1.2534363412216253633861")
# 2) Phi Watson check
n = mpf(10)**4
wat = sqrt(2)/n - sqrt(2)/(2*n**2) + 3*sqrt(2)/(4*n**3) - 7*sqrt(2)/(8*n**4)
print("Phi(1e4) - Watson(4 terms) =", nstr(Phi(n)-wat, 12), " expected O(n^-5) ~", nstr(mpf('2.5')/n**5, 6))
# product Φ/(1+e/2) n^{-4} coefficient = -3√2/2 : check numerically
prod = Phi(n)/(1+1/(2*n))
series = sqrt(2)/n - sqrt(2)/n**2 + 5*sqrt(2)/(4*n**3) - 3*sqrt(2)/(2*n**4)
print("Φ/(1+e/2) - series(4 terms) =", nstr(prod-series, 12))
# 3) A, B constants (prev step_07) at higher precision
def g(w):
    w = mpf(w); return (1 + w - sqrt(1+w*w))/w
A = quad(lambda w: g(w)-1, [0,1]) + quad(lambda w: g(w)-1/w, [1, inf])
B = -quad(lambda w: (g(w)-1)*log(w), [0,1]) - quad(lambda w: (g(w)-1/w)*log(w), [1, inf])
L2 = log(2); z2 = pi**2/6
print("A =", nstr(A, 40), "  A+log2 =", nstr(A+L2, 12))
print("B =", nstr(B, 40))
print("B - (z2/2 + L2^2/2 - L2) =", nstr(B - (z2/2+L2**2/2-L2), 12))
# 4) rho = 0.25 and rho=3 (via symmetry E_ρ(n)=E_{1/ρ}(ρn))
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
def R_rho_pred(n, rho):
    rho = mpf(rho); L = log(mpf(n))
    lam = log(1+rho**2)/2
    c1 = -(1+rho**3)/(rho*sqrt(1+rho**2))
    mt1 = 1-L2; mt2 = 1-L2+L2**2/2+z2/2
    return c1/rho - ((L+lam)**2/2 + mt1*(L+lam) + mt2)/rho
# rho=3 via E_{1/3}(3n)
for rho, n, nn in [('0.25', 4000, 4000), ('3', 1000, 3000)]:
    if mpf(rho) < 1:
        Ed = E_rho_direct(n, rho); Rnum = n*(n*Ed - sqrt(1+mpf(rho)**2)/mpf(rho)); Rp = R_rho_pred(n, rho)
    else:
        Ed = E_rho_direct(nn, 1/mpf(rho))  # E_{1/ρ}(ρn)
        Rnum = n*(n*Ed - sqrt(1+mpf(rho)**2)/mpf(rho)); Rp = R_rho_pred(n, rho)
    print("rho=%s n=%d: E=%s  Rnum=%s  Rpred=%s  diff=%s  (log^3 n)/n=%s" % (
        rho, n, nstr(Ed,25), nstr(Rnum,18), nstr(Rp,18), nstr(Rnum-Rp,8), nstr(log(n)**3/n,6)))
