# step_18.py : accurate C0, C' in sigma-form via GL; test poly and tail formulas
import math, numpy as np
from mpmath import mp, mpf, sqrt, exp, log, cosh, nstr, pi, zeta
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
            dP = deg*(p0 - x*p1)/(1-x*x); dx = p1/dP; x -= dx
            if abs(dx) < mpf(10)**(-mp.dps+10): break
        p0, p1 = mpf(1), x
        for k in range(2, deg+1):
            p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
        dP = deg*(p0 - x*p1)/(1-x*x)
        nodes.append((x, 2/((1-x*x)*dP*dP)))
    _gl_cache[deg] = nodes
    return nodes
def gli(func, a, b, deg=30):
    mid=(a+b)/2; half=(b-a)/2
    return sum(w*func(mid+half*x) for x,w in gl_nodes(deg))*half

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)

def C0_sig(n):  # ∫_0^∞ σ f(n√2 e^{-σ}) dσ
    n = mpf(n)
    sgmax = log(n*sqrt(2)) + 50
    npan = int(math.ceil(float(sgmax)/3.0))
    return sum(gli(lambda sg: sg*fW(n*sqrt(2)*exp(-sg)), 3*j, min(3*(j+1), sgmax)) for j in range(npan) if 3*j < sgmax)
def Cp_sig(n):  # ∫_0^∞ e^{-eσ} σ f(n√2 e^{-σ}) dσ
    n = mpf(n); e = 1/n
    sgmax = log(n*sqrt(2)) + 50
    npan = int(math.ceil(float(sgmax)/3.0))
    return sum(gli(lambda sg: exp(-e*sg)*sg*fW(n*sqrt(2)*exp(-sg)), 3*j, min(3*(j+1), sgmax)) for j in range(npan) if 3*j < sgmax)

L2 = log(2); z2 = pi**2/6; z3 = zeta(3); z4 = pi**4/90
e1=-L2; e2=L2**2/2+z2/2; e3=-z3/4-L2*z2/2-L2**3/6; e4=z4/4+L2*z3/4+L2**2*z2/4+z2**2/8+L2**4/24
eks=[mpf(1),e1,e2,e3,e4]; mt=[sum(eks[:k+1]) for k in range(5)]
def poly(k, N):
    L = log(mpf(N)); s = L**(k+1)/mpf(k+1)
    for i in range(1, k+2):
        s += mpf(math.factorial(k))/mpf(math.factorial(k+1-i))*mt[i]*L**(k+1-i)
    return s
for n in [1000, 1000000]:
    n_ = mpf(n); N = n_*sqrt(2)
    C0 = C0_sig(n_); Cp = Cp_sig(n_); D2 = Cp - C0
    p1 = poly(1, N); p2 = poly(2, N)
    print("n =", n)
    print("  C0(sig) =", nstr(C0, 45))
    print("  poly1   =", nstr(p1, 45))
    print("  C0-poly1 =", nstr(C0-p1, 15), "  vs -1/(2N) =", nstr(-1/(2*N), 15))
    print("  D2 = Cp-C0 =", nstr(D2, 25))
    print("  -e*poly2 =", nstr(-p2/n_, 25), "  D2 - (-e poly2) =", nstr(D2+p2/n_, 15))
    # poly series for D2
    D2poly = sum((-1/n_)**j/mpf(math.factorial(j))*poly(j+1, N) for j in range(1, 5))
    print("  D2poly(j=1..4) =", nstr(D2poly, 25), "  D2-D2poly =", nstr(D2-D2poly, 15))
