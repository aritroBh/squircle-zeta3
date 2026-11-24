"""
fnlib.py -- core library for Ahmed's Open Problem 3 (length of the graph of f^n).

Notation (see result.md):
  y_n = f^n on [0,1],  D_n = 2 - length(graph y_n),  mu = 1/(n-1),  a = f'(1).
  g(s)  = 1 + s - sqrt(1+s^2)            (taxicab defect density)
  g'(s) = 1 - s/sqrt(1+s^2)              (a probability density on (0,inf))
  R(w)  = Gamma(1+w/2) Gamma((1-w)/2) / (sqrt(pi) (1+w)) = int_0^inf g'(s) s^w ds
  A(t)  = -log f(1-t),   B(t) = log f'(1-t)
  T(w,mu): solution of A(T) - mu B(T) = w   (double power series)
  Expansion: D_n ~ sum_{i,j} t_{ij} mu^{i+j} [E_i(n) - Tail_i(n)]
     E_i(n)    = d^i/dz^i [ n^z R(-z) ]_{z=0}
     Tail_i(n) = int_{na}^inf g'(s) log(n/s)^i ds
"""
import mpmath as mp
import sympy as sp

# ----------------------------------------------------------------- basic functions
def g(s):
    """cancellation-free 1+s-sqrt(1+s^2)"""
    return 2*s/(1 + s + mp.sqrt(1 + s*s))

def gprime(s):
    r = mp.sqrt(1 + s*s)
    return 1/(r*(r + s))

def R(w):
    return mp.gamma(1 + w/2)*mp.gamma((1 - w)/2)/(mp.sqrt(mp.pi)*(1 + w))

def logR_coeffs(K):
    """Taylor coefficients a_j of log R(eps) = sum_j a_j eps^j, j=0..K (exact formula)."""
    c = [mp.mpf(0)]*(K + 1)
    if K >= 1:
        c[1] = mp.log(2) - 1
    for j in range(2, K + 1):
        if j % 2 == 0:
            c[j] = (mp.zeta(j) + 1)/j
        else:
            c[j] = (mp.zeta(j)*(1 - mp.mpf(2)**(1 - j)) - 1)/j
    return c

# ----------------------------------------------------------------- univariate power series (lists)
def ps_mul(p, q, K):
    r = [mp.mpf(0)]*(K + 1)
    for i, pi in enumerate(p[:K + 1]):
        if pi == 0:
            continue
        for j, qj in enumerate(q[:K + 1 - i]):
            r[i + j] += pi*qj
    return r

def ps_exp(p, K):
    """exp of a power series with p[0] arbitrary"""
    p = list(p) + [mp.mpf(0)]*(K + 1 - len(p))
    e0 = mp.exp(p[0])
    # r' = p' r
    r = [mp.mpf(0)]*(K + 1)
    r[0] = mp.mpf(1)
    for k in range(1, K + 1):
        s = mp.mpf(0)
        for j in range(1, k + 1):
            s += j*p[j]*r[k - j]
        r[k] = s/k
    return [e0*x for x in r]

def ps_log(p, K):
    """log of a power series with p[0] > 0"""
    p = list(p) + [mp.mpf(0)]*(K + 1 - len(p))
    p0 = p[0]
    q = [x/p0 for x in p]
    r = [mp.mpf(0)]*(K + 1)
    r[0] = mp.log(p0)
    # (log q)' = q'/q  ->  q * r' = q'
    for k in range(1, K + 1):
        s = k*q[k]
        for j in range(1, k):
            s -= j*r[j]*q[k - j]
        r[k] = s/k
    return r

def moments_xi(K):
    """rho_m = E[xi^m], xi = log(1/S), S ~ g'(s) ds;  generating function R(-z) = sum rho_m z^m/m!"""
    a = logR_coeffs(K)
    am = [a[j]*(-1)**j for j in range(K + 1)]   # log R(-z)
    Rm = ps_exp(am, K)                            # R(-z) Taylor coefficients
    return [Rm[m]*mp.factorial(m) for m in range(K + 1)]

# ----------------------------------------------------------------- bivariate series T(w,mu)
# store as dict (i,j) -> coeff for w^i mu^j, total degree <= K
def bs_mul(P, Q, K):
    out = {}
    for (i1, j1), c1 in P.items():
        for (i2, j2), c2 in Q.items():
            if i1 + j1 + i2 + j2 <= K:
                key = (i1 + i2, j1 + j2)
                out[key] = out.get(key, 0) + c1*c2
    return out

def bs_compose(coeffs, T, K):
    """sum_k coeffs[k] T^k, T without constant term"""
    out = {}
    Tp = {(0, 0): mp.mpf(1)}
    for k, ck in enumerate(coeffs[:K + 1]):
        if k > 0:
            Tp = bs_mul(Tp, T, K)
        if ck != 0:
            for key, v in Tp.items():
                out[key] = out.get(key, 0) + ck*v
    return out

def solve_T(Acoef, Bcoef, K):
    """A(T) - mu*B(T) = w,  A(0)=0, A'(0)=alpha1 != 0."""
    a1 = Acoef[1]
    T = {}
    for it in range(K + 2):
        AT = bs_compose([0, 0] + list(Acoef[2:]), T, K)    # A(T) - a1 T
        BT = bs_compose(Bcoef, T, K)
        new = {(1, 0): mp.mpf(1)/a1}
        for (i, j), v in BT.items():
            if i + j + 1 <= K:
                new[(i, j + 1)] = new.get((i, j + 1), 0) + v/a1
        for key, v in AT.items():
            new[key] = new.get(key, 0) - v/a1
        T = {k: v for k, v in new.items() if k[0] + k[1] <= K}
    return T

# ----------------------------------------------------------------- expansion polynomials P_k(L)
def tail_coeffs(i, a, kmax):
    """Tail_i(n) = sum_{k>=1} c_{i,k} (n a)^{1-2k}  (exact for n a > 1)"""
    la = mp.log(a)
    out = []
    for k in range(1, kmax + 1):
        bk = -mp.binomial(mp.mpf(-1)/2, k)       # g'(s) = sum_k bk s^{-2k}
        sm = mp.mpf(0)
        for l in range(i + 1):
            sm += mp.binomial(i, l)*la**(i - l)*mp.factorial(l)/mp.mpf(2*k - 1)**(l + 1)
        out.append(bk*(-1)**i*sm)                # times (n a)^{1-2k}
    return out

def expansion_polys(Acoef, Bcoef, K):
    """Return list P[k] (k=0..K), each a list of coefficients in L=log n (degree <= k),
    such that D_n ~ sum_k mu^k P_k(log n)."""
    a = Acoef[1]           # = f'(1)
    T = solve_T(Acoef, Bcoef, K)
    rho = moments_xi(K)
    P = [[mp.mpf(0)]*(k + 1) for k in range(K + 1)]
    # main part: t_ij mu^{i+j} E_i, E_i = sum_m C(i,m) L^{i-m} rho_m
    for (i, j), t in T.items():
        k = i + j
        for m in range(i + 1):
            P[k][i - m] += t*mp.binomial(i, m)*rho[m]
    # tail part: - t_ij mu^{i+j} Tail_i(n);  (n a)^{1-2k} = a^{1-2k} mu^{2k-1} (1+mu)^{1-2k}
    for (i, j), t in T.items():
        k0 = i + j
        if k0 >= K:
            continue
        kmax = (K - k0 + 1)//2 + 1
        tc = tail_coeffs(i, a, kmax)
        for kk, c in enumerate(tc, start=1):
            p = 2*kk - 1
            # mu^{k0+p} (1+mu)^{-p} = sum_r C(-p,r) mu^{k0+p+r}
            for r in range(0, K - k0 - p + 1):
                P[k0 + p + r][0] -= t*c*a**(1 - 2*kk)*mp.binomial(-p, r)
    return P

def eval_expansion(P, n, K):
    n = mp.mpf(n)
    mu = 1/(n - 1)
    L = mp.log(n)
    s = mp.mpf(0)
    for k in range(1, K + 1):
        s += mu**k*mp.polyval(P[k][::-1], L)
    return s

# ----------------------------------------------------------------- test functions
x = sp.Symbol('x', positive=True)
TESTS = {
    'x':        x,
    'x^2':      x**2,
    'expm1':    (sp.exp(x) - 1)/(sp.E - 1),
    'x+x^3':    (x + x**3)/2,
    'mobius':   2*x/(1 + x),
    'sin4':     sp.sin(sp.pi*x/4)/sp.sin(sp.pi/4),
}

def jets(expr, K):
    """A(t) = -log f(1-t), B(t) = log f'(1-t) Taylor coefficients (numeric, current mp.dps)."""
    t = sp.Symbol('t')
    fexpr = expr.subs(x, 1 - t)
    fpexpr = sp.diff(expr, x).subs(x, 1 - t)
    fs = sp.series(fexpr, t, 0, K + 2).removeO()
    fps = sp.series(fpexpr, t, 0, K + 2).removeO()
    fc = [mp.mpf(sp.N(fs.coeff(t, k), mp.mp.dps + 20)) for k in range(K + 2)]
    fpc = [mp.mpf(sp.N(fps.coeff(t, k), mp.mp.dps + 20)) for k in range(K + 2)]
    A = [-c for c in ps_log(fc, K + 1)]
    B = ps_log(fpc, K + 1)
    return A, B

def make_callables(expr):
    f = sp.lambdify(x, expr, 'mpmath')
    fp = sp.lambdify(x, sp.diff(expr, x), 'mpmath')
    return f, fp

# ----------------------------------------------------------------- high precision defect D_n
def defect(f, fp, n, a, dps_target=50, vtail=None, verbose=False):
    """D_n = int_0^1 g(n f^{n-1} f') dx computed in tau = 1-x = v/n.
    Returns (D, rigorous bound on the neglected far tail)."""
    n = mp.mpf(n)
    if vtail is None:
        vtail = (dps_target + 15)*mp.log(10)/a + 2*mp.log(n*a)/a + 50
    vmax = min(vtail, n)              # tau <= 1
    nm1 = n - 1

    def integrand(v):
        xx = 1 - v/n
        fx = f(xx)
        if fx <= 0:
            return mp.mpf(0)
        s = n*mp.exp(nm1*mp.log(fx))*fp(xx)
        return g(s)
    c = mp.log(n*a)/a
    pts = set([mp.mpf(0), vmax])
    for k in [-8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
        p = c + mp.mpf(k)/a
        if 0 < p < vmax:
            pts.add(p)
    pts = sorted(pts)
    val = mp.quad(integrand, pts)/n
    # rigorous tail bound: int_{x< 1 - vmax/n} g(y') dx <= y(1 - vmax/n) = f^n
    if vmax < n:
        bound = f(1 - vmax/n)**n
    else:
        bound = mp.mpf(0)
    return val, bound
