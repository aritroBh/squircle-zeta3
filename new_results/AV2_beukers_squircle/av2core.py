
"""Core functions for AV2: Beukers zeta(3) integrals on squircle/superball domains."""
from mpmath import (mp, mpf, loggamma, digamma, exp, log, log1p, polylog, pi,
                    zeta, quad, fabs, inf, nsum, sin, cos, tan, sqrt, euler)

def phi(w):
    """phi(w) = w^3 Gamma(w)^2/Gamma(2w+1) * (psi(2w+1)-psi(w));  phi(0)=1.
    B2(p) = sum_{n>=1} phi(n*eps)/n^3 with eps=1/p."""
    w = mpf(w)
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w))

def B2_term(n, eps):
    w = n*eps
    return w**3 * exp(2*loggamma(w) - loggamma(2*w+1)) * (digamma(2*w+1) - digamma(w)) / mpf(n)**3

def B2_series(p, tol=None, nmin_mult=30):
    """Exact Dirichlet series: B2(p) = sum phi(n/p)/n^3. Geometric tail ~4^(-n/p)."""
    p = mpf(p); eps = 1/p
    if tol is None:
        tol = mpf(10)**(-(mp.dps - 10))
    S = mpf(0); n = 1
    nmin = int(nmin_mult*p) + 10
    while True:
        T = B2_term(n, eps)
        S += T
        if n > nmin and fabs(T) < tol:
            break
        n += 1
    return S, n

def B2_quad(p):
    """1D integral: B2(p) = (1/p) int_0^{pi/4} (cos sin)^{2/p-1} g(c) dth, c=(cos sin)^{2/p},
    g(c) = (Li_2(c) + ln c ln(1-c))/c = (pi^2/6 - Li_2(1-c))/c.   t = tan(theta)."""
    p = mpf(p); eps = 1/p; a = 2*eps
    def f(t):
        if t == 0:
            return mpf(0)
        u = t/(1+t*t)
        c = u**a
        g = (polylog(2, c) + log(c)*log1p(-c))/c
        return u**(a-1) * g / (1+t*t)
    return eps*quad(f, [0, 1])

def I_moment(a, b, p):
    """Dirichlet moment: int_{D_p} x^a y^b dxdy = (1/p^2) G((a+1)/p) G((b+1)/p)/G((a+b+2)/p+1)."""
    p = mpf(p)
    return exp(2*log(1/p) + loggamma((a+1)/p) + loggamma((b+1)/p) - loggamma((a+b+2)/p + 1))

def M_moment(j, k, p, tol=None, mmin=None):
    """M(j,k,p) = int_{D_p} -log(xy)/(1-xy) x^j y^k dxdy
    = -(1/p) sum_{m>=0} I(j+m,k+m,p) [psi((j+m+1)/p)+psi((k+m+1)/p)-2 psi((j+k+2m+2)/p+1)]."""
    p = mpf(p); ep = 1/p
    if tol is None:
        tol = mpf(10)**(-(mp.dps - 10))
    S = mpf(0); m = 0
    if mmin is None:
        mmin = int(30*p) + 10
    while True:
        a1 = (j+m+1)*ep; b1 = (k+m+1)*ep; c1 = (j+k+2*m+2)*ep + 1
        Im = exp(2*log(ep) + loggamma(a1) + loggamma(b1) - loggamma(c1))
        T = -ep * Im * (digamma(a1) + digamma(b1) - 2*digamma(c1))
        S += T
        if m > mmin and fabs(T) < tol:
            break
        m += 1
    return S

def shifted_legendre_coeffs(n):
    """Shifted Legendre P_n(x) = sum_k c_k x^k on [0,1], P_n(x)=(1/n!)(d/dx)^n x^n(1-x)^n."""
    from math import comb
    return [(-1)**(n-k) * comb(n, k) * comb(n+k, k) for k in range(n+1)]

def L_legendre(n, p, Mdict=None):
    """L_n(p) = int_{D_p} -log(xy)/(1-xy) P_n(x) P_n(y) dxdy via moment sums."""
    cs = shifted_legendre_coeffs(n)
    tot = mpf(0)
    for i, ci in enumerate(cs):
        for j, cj in enumerate(cs):
            if ci == 0 or cj == 0:
                continue
            key = (min(i, j), max(i, j))
            if Mdict is not None and key in Mdict:
                M = Mdict[key]
            else:
                M = M_moment(i, j, p)
                if Mdict is not None:
                    Mdict[key] = M
            tot += ci*cj*M
    return tot



def B2_quad2(p):
    """1D form B2(p) = (1/2) int_0^inf u^{2/p} g(c) dv with u=1/(2 cosh(pv/2)), c=u^{2/p},
    g(c)=(Li_2(c)+ln c ln(1-c))/c.  Integrand ~ e^{-v}(1+v) at infinity: smooth for all p."""
    p = mpf(p); ep = 1/p
    half = mpf(1)/2
    two_ep = 2*ep
    def f(v):
        u = 1/(2*mp.cosh(p*v*half))
        c = u**two_ep
        g = (polylog(2, c) + log(c)*log1p(-c))/c
        return c * g      # u^{2/p} = c
    V = (mp.dps + 25)*log(10)   # tail ~ e^{-V}(1+V)
    return half*quad(f, [0, V])

def B2_quad3(p):
    """1D form B2(p) = (1/2) int_0^{2^{-2/p}} g(c)/sqrt(1-4 c^p) dc  (third method)."""
    p = mpf(p); ep = 1/p
    cup = 2**(-2*ep)
    def h(c):
        if c == 0:
            return mpf(0)
        g = (polylog(2, c) + log(c)*log1p(-c))/c
        return g/mp.sqrt(1-4*c**p)
    return mpf(1)/2*quad(h, [0, cup])

def C2_direct(W=None):
    """C2 = int_0^inf (phi(w)-1) w^{-3} dw = int_0^W (...)dw - 1/(2W^2) + O(4^{-W}).
    Near w=0 use Taylor phi(w)-1 = -zeta(3) w^3 + (9/4) zeta(4) w^4 + ..."""
    if W is None:
        W = mpf(200)
    z3 = zeta(3); c4 = mpf(9)/4*zeta(4); c5 = 3*zeta(2)*z3 - 9*zeta(5)
    def h(w):
        if w < mpf('1e-20'):
            return z3*0 - z3 + c4*w - c5*0 + c5*w**2  # -z3 + phi4 w + phi5 w^2
        return (phi(w)-1)/w**3
    val = quad(h, [0, 1]) + quad(h, [1, 20]) + quad(h, [20, W])
    return val - 1/(2*W*W)


def B3_quadA(p):
    """2D: x=s, y=t(1-s^p)^{1/p}; corner-stable evaluation:
    1-(1-xy)Z = [(1-A)+A(1-B)] + xy*Z  (all terms >= 0, no cancellation),
    A=(1-s^p)^{1/p}, B=(1-t^p)^{1/p}, Z=AB, xy=s*t*A;  q=1-xy >= 1-2^{-2/p} > 0."""
    p = mpf(p); ep = 1/p
    def f(s, t):
        lA = ep*log1p(-s**p)
        lB = ep*log1p(-t**p)
        A = exp(lA); B = exp(lB)
        one_A = -mp.expm1(lA); one_B = -mp.expm1(lB)
        xy = s*t*A
        val = (one_A + A*one_B) + xy*(A*B)
        return A * (-log(val))/(1-xy)
    return quad(f, [0, 1], [0, 1])

def B3_quadB(p):
    """2D p-polar: corner-stable: 1-(1-r^2 c)Q = (1-Q) + r^2 c Q."""
    p = mpf(p); ep = 1/p; a = 2*ep
    def f(th, r):
        cs = cos(th)*sin(th)
        c = cs**a
        lQ = ep*log1p(-r**p)
        Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        q = 1-r*r*c
        return (2*ep) * cs**(a-1) * r * (-log(val))/q
    return quad(f, [0, pi/2], [0, 1])


# ---------------- custom tensor tanh-sinh quadrature ----------------
def _ts_nodes(a, b, h, n):
    """nodes/weights for int_a^b f dx via tanh-sinh: x = mid + half*tanh(pi/2 sinh(kh))."""
    a = mpf(str(a)); b = mpf(str(b)); h = mpf(str(h))
    mid = (a+b)/2; half = (b-a)/2
    pih = pi*h/2
    out = []
    for k in range(-n, n+1):
        t = k*h
        sh = mp.sinh(t)
        u = mp.tanh(pih*0 + (pi/2)*sh) if False else mp.tanh((pi/2)*sh)
        w = (pi/2)*mp.cosh(t)/(mp.cosh((pi/2)*sh))**2
        out.append((mid + half*u, half*w*h))
    return out

def ts1(f, a, b, h, n):
    s = mpf(0)
    for x, w in _ts_nodes(a, b, h, n):
        s += w*f(x)
    return s

def ts2(f, a, b, c, d, h, n):
    xs = _ts_nodes(a, b, h, n)
    ys = _ts_nodes(c, d, h, n)
    s = mpf(0)
    for x, wx in xs:
        for y, wy in ys:
            s += wx*wy*f(x, y)
    return s

def B3_tsA(p, h, n):
    """square-mapped 2D with stable corner."""
    p = mpf(p); ep = 1/p
    def f(s, t):
        lA = ep*log1p(-s**p); lB = ep*log1p(-t**p)
        A = exp(lA); B = exp(lB)
        one_A = -mp.expm1(lA); one_B = -mp.expm1(lB)
        xy = s*t*A
        val = (one_A + A*one_B) + xy*(A*B)
        if val <= 0:
            return mpf(0)   # underflowed corner node; its weight is ~0 anyway
        return A * (-log(val))/(1-xy)
    return ts2(f, 0, 1, 0, 1, h, n)

def B3_tsB(p, h, n):
    """p-polar 2D with stable corner; th in [0, pi/2], r in [0,1]."""
    p = mpf(p); ep = 1/p; a = 2*ep
    hp = pi/2
    def f(th, r):
        cs = cos(th)*sin(th)
        if cs < 0:
            cs = mpf(0)   # endpoint rounding artifact
        c = cs**a
        lQ = ep*log1p(-r**p)
        Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        if val <= 0:
            return mpf(0)
        q = 1-r*r*c
        return (2*ep) * cs**(a-1) * r * (-log(val))/q
    return ts2(f, 0, hp, 0, 1, h, n)


def F_cp(c, p, h2, n2):
    """F(c,p) = int_0^1 r (-log(1-(1-r^2 c)(1-r^p)^{1/p}))/(1-r^2 c) dr.  F(c,inf)=g(c)/2."""
    p = mpf(p); ep = 1/p
    def f(r):
        lQ = ep*log1p(-r**p)
        Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        if val <= 0:
            return mpf(0)
        q = 1-r*r*c
        return r*(-log(val))/q
    return ts1(f, 0, 1, h2, n2)

def B3_cform(p, h1, n1, h2, n2):
    """B3(p) = 2 int_0^{2^{-2/p}} F(c,p)/sqrt(1-4 c^p) dc   (outer 1D TS x inner 1D TS)."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    def outer(c):
        if c == 0:
            return mpf(0)
        den = 1-4*c**p
        if den <= 0:
            return mpf(0)   # endpoint node; weight ~ 0
        return F_cp(c, p, h2, n2)/mp.sqrt(den)
    return 2*ts1(outer, 0, cstar, h1, n1)

def B3_vform(p, h1, n1, h2, n2):
    """B3(p) = 2 int_0^inf c F(c,p) dv,  c = (2 cosh(pv/2))^{-2/p}."""
    p = mpf(p); ep = 1/p
    half = mpf(1)/2
    def outer(v):
        u = 1/(2*mp.cosh(p*v*half))
        c = u**(2*ep)
        return c * F_cp(c, p, h2, n2)
    V = (mp.dps + 40)*log(10)
    tot = mpf(0)
    for (a, b) in [(0, 16), (16, 64), (64, V)]:
        tot += ts1(outer, a, b, h1, n1)
    return 2*tot


def F_cp2(c, p, h2, n2):
    """F(c,p) via r=(1-e^{-s})^{1/p}: F = (1/p) int_0^inf (1-e^{-s})^{2/p-1} e^{-s} (-log val)/q ds,
    val = (1-e^{-s/p}) + e^{-s/p} (1-e^{-s})^{2/p} c  (all terms >= 0),  q = 1-(1-e^{-s})^{2/p} c."""
    p = mpf(p); ep = 1/p
    def f(s):
        if s == 0:
            return mpf(0)
        ems = -mp.expm1(-s)          # 1-e^{-s}
        r2 = ems**(2*ep)             # r^2 = (1-e^{-s})^{2/p}
        q = 1-r2*c
        one_Q = -mp.expm1(-s*ep)     # 1-e^{-s/p}
        val = one_Q + exp(-s*ep)*(ems**(2*ep))*c
        if val <= 0:
            return mpf(0)
        return ep * ems**(2*ep-1) * exp(-s) * (-log(val))/q
    S = (mp.dps + 25)*log(10)
    tot = mpf(0)
    for (a, b) in [(0, 12), (12, 48), (48, S)]:
        tot += ts1(f, a, b, h2, n2)
    return tot

def B3_cform2(p, h1, n1, h2, n2):
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    def outer(c):
        if c == 0:
            return mpf(0)
        den = 1-4*c**p
        if den <= 0:
            return mpf(0)
        return F_cp2(c, p, h2, n2)/mp.sqrt(den)
    return 2*ts1(outer, 0, cstar, h1, n1)


def F_cp3(c, p, h2, n2):
    """F(c,p) = int_0^1 r (-log val)/q dr, split at r=1/2;
    piece 2 uses r = 1-rho/p, rho in [0, p/2]: resolves the (1-r)^{1/p} layer."""
    p = mpf(p); ep = 1/p
    R = mpf(1)/2
    def f1(r):
        lQ = ep*log1p(-r**p)
        Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        if val <= 0:
            return mpf(0)
        return r*(-log(val))/(1-r*r*c)
    L = p/2   # rho in [0, p/2] <=> r in [1/2, 1]
    def f2(rho):
        r = 1-rho/p
        om = -mp.expm1(p*log1p(-rho/p))   # 1-(1-rho/p)^p
        Q = om**ep
        val = (1-Q) + (r*r*c)*Q
        if val <= 0:
            return mpf(0)
        return (1/p)*r*(-log(val))/(1-r*r*c)
    return ts1(f1, 0, R, h2, n2) + ts1(f2, 0, L, h2, n2)

def B3_cform3(p, h1, n1, h2, n2):
    """B3(p) = 2 int_0^{c*} F(c,p)/sqrt(1-4c^p) dc with split inner integral."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    def outer(c):
        if c == 0:
            return mpf(0)
        den = 1-4*c**p
        if den <= 0:
            return mpf(0)
        return F_cp3(c, p, h2, n2)/mp.sqrt(den)
    return 2*ts1(outer, 0, cstar, h1, n1)

def B3_vform3(p, h1, n1, h2, n2):
    """B3(p) = 2 int_0^inf c F(c,p) dv with split inner integral."""
    p = mpf(p); ep = 1/p
    half = mpf(1)/2
    def outer(v):
        u = 1/(2*mp.cosh(p*v*half))
        c = u**(2*ep)
        return c * F_cp3(c, p, h2, n2)
    V = (mp.dps + 40)*log(10)
    tot = mpf(0)
    for (a, b) in [(0, 16), (16, 64), (64, V)]:
        tot += ts1(outer, a, b, h1, n1)
    return 2*tot


def B3_cform4(p, h1, n1, h2, n2, G=30):
    """Outer split: [0, c*e^{-G/p}] plain + c=c*e^{-gam/p}, gam in [0,G]."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    c0 = cstar*exp(-G*ep)
    def F(c):
        return F_cp3(c, p, h2, n2)
    def out1(c):
        if c == 0:
            return mpf(0)
        return F(c)/mp.sqrt(1-4*c**p)
    def out2(gam):
        c = cstar*exp(-gam*ep)
        den = -mp.expm1(-gam)     # 1-e^{-gam}
        if den <= 0:
            return mpf(0)
        return F(c)/mp.sqrt(den) * (c*ep)
    return 2*(ts1(out1, 0, c0, h1, n1) + ts1(out2, 0, G, h1, n1))


def F_cpQ(c, p, h2, n2):
    """F(c,p) = int_0^1 Q^{p-1} (1-Q^p)^{2/p-1} (-log val)/q dQ,
    val = (1-Q) + Q (1-Q^p)^{2/p} c  (>=0, no cancellation),  q = 1-(1-Q^p)^{2/p} c."""
    p = mpf(p); ep = 1/p
    def f(Q):
        om = 1-Q**p if Q < 1 else mpf(0)     # 1-Q^p
        r2 = om**(2*ep)                      # r^2 = (1-Q^p)^{2/p}
        q = 1-r2*c
        val = (1-Q) + Q*r2*c
        if val <= 0:
            return mpf(0)
        return Q**(p-1) * om**(2*ep-1) * (-log(val))/q
    return ts1(f, 0, 1, h2, n2)

def B3_cformQ(p, h1, n1, h2, n2):
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    def outer(c):
        if c == 0:
            return mpf(0)
        den = 1-4*c**p
        if den <= 0:
            return mpf(0)
        return F_cpQ(c, p, h2, n2)/mp.sqrt(den)
    return 2*ts1(outer, 0, cstar, h1, n1)

def B3_vformQ(p, h1, n1, h2, n2):
    p = mpf(p); ep = 1/p
    half = mpf(1)/2
    def outer(v):
        u = 1/(2*mp.cosh(p*v*half))
        c = u**(2*ep)
        return c * F_cpQ(c, p, h2, n2)
    V = (mp.dps + 40)*log(10)
    tot = mpf(0)
    for (a, b) in [(0, 16), (16, 64), (64, V)]:
        tot += ts1(outer, a, b, h1, n1)
    return 2*tot


def F_cpU(c, p, h2, n2):
    """F(c,p) = int_0^inf e^{-2u} (-log val)/q du,  r=e^{-u},
    val = (1-Q) + e^{-2u} c Q  (>=0),  Q=(1-e^{-pu})^{1/p},  q=1-e^{-2u}c."""
    p = mpf(p); ep = 1/p
    def f(u):
        if u == 0:
            return mpf(0)
        om = -mp.expm1(-p*u)          # 1-e^{-pu}
        Q = om**ep
        val = (1-Q) + exp(-2*u)*c*Q
        if val <= 0:
            return mpf(0)
        return exp(-2*u)*(-log(val))/(1-exp(-2*u)*c)
    S = ((mp.dps + 30)*log(10))/2
    tot = mpf(0)
    for (a, b) in [(0, 10), (10, 40), (40, S)]:
        tot += ts1(f, a, b, h2, n2)
    return tot

def B3_cformU(p, h1, n1, h2, n2):
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    def outer(c):
        if c == 0:
            return mpf(0)
        den = 1-4*c**p
        if den <= 0:
            return mpf(0)
        return F_cpU(c, p, h2, n2)/mp.sqrt(den)
    return 2*ts1(outer, 0, cstar, h1, n1)


def F_cpOM(c, p, h2, n2):
    """F(c,p) = (1/p) int_0^1 (1-om)^{2/p-1} (-log val)/q dom,  om=1-r^p,
    val = (1-Q) + Q*(1-om)^{2/p} c with Q = om^{1/p}  (>=0),  q = 1-(1-om)^{2/p} c."""
    p = mpf(p); ep = 1/p
    def f(om):
        if om == 0:
            return mpf(0)
        Q = om**ep
        r2 = (1-om)**(2*ep)
        q = 1-r2*c
        val = (1-Q) + Q*r2*c
        if val <= 0:
            return mpf(0)
        return ep * (1-om)**(2*ep-1) * (-log(val))/q
    return ts1(f, 0, 1, h2, n2)

def F_cp_split(c, p, hA, nA, hB, nB):
    """F = int_0^{0.9} f(r) dr  +  (1/p) int_0^{1-0.9^p} (1-om)^{2/p-1}(-log val)/q dom."""
    p = mpf(p); ep = 1/p
    def f1(r):
        lQ = ep*log1p(-r**p); Q = exp(lQ)
        one_Q = -mp.expm1(lQ)
        val = one_Q + (r*r*c)*Q
        if val <= 0: return mpf(0)
        return r*(-log(val))/(1-r*r*c)
    om1 = 1-mpf('0.9')**p
    def f2(om):
        if om == 0: return mpf(0)
        Q = om**ep
        r2 = (1-om)**(2*ep)
        q = 1-r2*c
        val = (1-Q) + Q*r2*c
        if val <= 0: return mpf(0)
        return ep * (1-om)**(2*ep-1) * (-log(val))/q
    return ts1(f1, 0, mpf('0.9'), hA, nA) + ts1(f2, 0, om1, hB, nB)


class FPrecomp:
    """Fast F(c,p) evaluation: precompute r-grid with Q factors (c-independent).
    F(c) = sum_k w_k r_k (-log(oneQ_k + r2_k c Q_k))/(1-r2_k c)."""
    def __init__(self, p, h2, n2):
        p = mpf(p); ep = 1/p
        self.grid = []
        for (r, w) in _ts_nodes(0, 1, h2, n2):
            lQ = ep*log1p(-r**p)
            Q = exp(lQ)
            oneQ = -mp.expm1(lQ)
            self.grid.append((r, w, r*r, Q, oneQ))
    def F(self, c):
        s = mpf(0)
        for (r, w, r2, Q, oneQ) in self.grid:
            val = oneQ + r2*c*Q
            if val <= 0:
                continue   # saturated node; true contribution is w*0*log = 0
            s += w*r*(-log(val))/(1-r2*c)
        return s

def B3_fast(p, h1, n1, FP, G=30):
    """Outer: [0, c*e^{-G/p}] + gamma-piece, using precomputed inner F."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    c0 = cstar*exp(-G*ep)
    def out1(c):
        if c == 0:
            return mpf(0)
        return FP.F(c)/mp.sqrt(1-4*c**p)
    def out2(gam):
        c = cstar*exp(-gam*ep)
        den = -mp.expm1(-gam)
        if den <= 0:
            return mpf(0)
        return FP.F(c)/mp.sqrt(den)*(c*ep)
    return 2*(ts1(out1, 0, c0, h1, n1) + ts1(out2, 0, G, h1, n1))


def B3_fast2(p, h1, n1, FP, G=40):
    """B3 = 2 F(c*) c* (1/p) B(1/p,1/2)  +  2 int [F(c)-F(c*)]/sqrt(1-4c^p) dc,
    second integral split at c0 = c* e^{-G/p} with gamma-piece c = c* e^{-gam/p}."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    Fstar = FP.F(cstar)
    beta = ep*exp(loggamma(ep) + loggamma(mpf(1)/2) - loggamma(ep + mpf(1)/2))
    base = 2*Fstar*cstar*beta
    c0 = cstar*exp(-G*ep)
    def out1(c):
        if c == 0:
            return mpf(0)
        return (FP.F(c)-Fstar)/mp.sqrt(1-4*c**p)
    def out2(gam):
        c = cstar*exp(-gam*ep)
        den = -mp.expm1(-gam)
        if den <= 0:
            return mpf(0)
        return (FP.F(c)-Fstar)/mp.sqrt(den)*(c*ep)
    return base + 2*(ts1(out1, 0, c0, h1, n1) + ts1(out2, 0, G, h1, n1))


class FPrecompC:
    """F(c,p) with cache: reuses evaluations across nested outer levels."""
    def __init__(self, p, h2, n2):
        p = mpf(p); ep = 1/p
        self.grid = []
        for (r, w) in _ts_nodes(0, 1, h2, n2):
            lQ = ep*log1p(-r**p)
            Q = exp(lQ)
            oneQ = -mp.expm1(lQ)
            self.grid.append((r, w, r*r, Q, oneQ))
        self.cache = {}
    def F(self, c):
        c = +c
        if c in self.cache:
            return self.cache[c]
        s = mpf(0)
        for (r, w, r2, Q, oneQ) in self.grid:
            val = oneQ + r2*c*Q
            if val <= 0:
                continue
            s += w*r*(-log(val))/(1-r2*c)
        self.cache[c] = s
        return s

def B3_prod(p, h1, n1, FP, G=40):
    """cform with peak subtraction; levels nested (call twice with h1, h1/2 for reuse)."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    Fstar = FP.F(cstar)
    beta = ep*exp(loggamma(ep) + loggamma(mpf(1)/2) - loggamma(ep + mpf(1)/2))
    base = 2*Fstar*cstar*beta
    c0 = cstar*exp(-G*ep)
    def out1(c):
        if c == 0:
            return mpf(0)
        return (FP.F(c)-Fstar)/mp.sqrt(1-4*c**p)
    def out2(gam):
        c = cstar*exp(-gam*ep)
        den = -mp.expm1(-gam)
        if den <= 0:
            return mpf(0)
        return (FP.F(c)-Fstar)/mp.sqrt(den)*(c*ep)
    return base + 2*(ts1(out1, 0, c0, h1, n1) + ts1(out2, 0, G, h1, n1)), base, Fstar

def B3_vprod(p, h1, n1, FP, base=None, Fstar=None):
    """vform with peak subtraction: B3 = base + 2 int_0^inf c (F(c)-F(c*)) dv."""
    p = mpf(p); ep = 1/p
    cstar = 2**(-2*ep)
    if Fstar is None:
        Fstar = FP.F(cstar)
    if base is None:
        beta = ep*exp(loggamma(ep) + loggamma(mpf(1)/2) - loggamma(ep + mpf(1)/2))
        base = 2*Fstar*cstar*beta
    half = mpf(1)/2
    def outer(v):
        u = 1/(2*mp.cosh(p*v*half))
        c = u**(2*ep)
        return c*(FP.F(c)-Fstar)
    V = (mp.dps + 40)*log(10)
    tot = mpf(0)
    for (a, b) in [(0, 16), (16, 64), (64, V)]:
        tot += ts1(outer, a, b, h1, n1)
    return base + 2*tot
