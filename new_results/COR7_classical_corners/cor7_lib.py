"""cor7_lib.py -- shared machinery for COR7.
M3 formula (G1, B3): M3 = -zeta(3) Phi(c) + zeta(2) Psi(c) + R3(c), proved to all orders fitted;
weight split as in COR5: W4 = R3^(4) + z2 Psi^(2) + z3 * ((1-c^2)/(2c)) atanh c,
                         W3 = R3^(3) + z2 Psi^(1) - z3/2,
                         W2 = R3^(2) + z2/2.
HPLs evaluated by the exact recursive series (COR4 hpl_series.py).
"""
import json, math, re
import mpmath as mp

BASE = '<project>/research/'
G1 = BASE + '10_frontier/G1_zeta3_any_angle/'

# ---------------- HPL series (verbatim from COR4 hpl_series.py) ----------------
def hpl_coeffs(w, N):
    if len(w) == 0:
        c = [mp.mpf(0)] * (N + 1); c[0] = mp.mpf(1); return c
    a, rest = w[0], w[1:]
    c = hpl_coeffs(rest, N)
    out = [mp.mpf(0)] * (N + 1)
    if a == 0:
        for n in range(1, N + 1):
            out[n] = c[n] / n
    elif a == 1:
        S = mp.mpf(0)
        for n in range(1, N + 1):
            S += c[n - 1]
            out[n] = S / n
    else:
        A = mp.mpf(0)
        for n in range(1, N + 1):
            A = c[n - 1] - A
            out[n] = A / n
    return out

def H_series(w, x, N=None):
    w = tuple(w)
    if w == (0,):
        return mp.log(x)
    if w == (1,):
        return -mp.log(1 - x)
    if w == (-1,):
        return mp.log(1 + x)
    if w == (0, 1):
        return mp.polylog(2, x)
    if w == (0, -1):
        return -mp.polylog(2, -x)
    if w == (1, 1):
        return mp.log(1 - x)**2/2
    if w == (-1, -1):
        return mp.log(1 + x)**2/2
    if w == (0, 0, 1):
        return mp.polylog(3, x)
    if w == (0, 0, -1):
        return -mp.polylog(3, -x)
    if w == (0, 0, 0, 1):
        return mp.polylog(4, x)
    if w == (0, 0, 0, -1):
        return -mp.polylog(4, -x)
    if N is None:
        N = int((mp.mp.dps + 15) / abs(math.log10(abs(x)))) + 10
    c = hpl_coeffs(w, N)
    tot = mp.mpf(0); xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn; xn *= x
    return tot

# ---------------- the fitted M3 data ----------------
def _parse(name, coef):
    num, den = (coef.split("/") + ["1"])[:2]
    v = mp.mpf(int(num))/int(den)
    if name == '1':
        return 0, (), v
    if name == 'c':
        return 1, (), v
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2) else ()
    return e, w, v

PSI = [_parse(n, c) for n, c in json.load(open(G1 + 'hpl_fit_k3_(1,0,0).json'))]
R3T = [_parse(n, c) for n, c in json.load(open(G1 + 'hpl_fit_k3_(0,0,0).json'))]

def _combo(terms, c, wts=None):
    tot = mp.mpf(0)
    for e, w, v in terms:
        if wts is not None and len(w) not in wts:
            continue
        tot += v * c**e * (mp.mpf(1) if not w else H_series(w, c))
    return tot

def Phi(c):
    return mp.mpf(1)/2 - (1 - c*c)/(2*c) * mp.atanh(c)

def M3_weights(c):
    """returns (M3, W4, W3, W2) at the current precision."""
    z2, z3 = mp.zeta(2), mp.zeta(3)
    R32 = _combo(R3T, c, wts={2}); R33 = _combo(R3T, c, wts={3}); R34 = _combo(R3T, c, wts={4})
    P0 = _combo(PSI, c, wts={0}); P1 = _combo(PSI, c, wts={1}); P2 = _combo(PSI, c, wts={2})
    W4 = R34 + z2*P2 + z3*(1 - c*c)/(2*c)*mp.atanh(c)
    W3 = R33 + z2*P1 - z3/2
    W2 = R32 + z2*P0
    return W4 + W3 + W2, W4, W3, W2

# ---------------- LLL relation finder (COR5 style) ----------------
def find_relation_rows(vals, dps, margin=30, delta=0.99):
    from flint import fmpz_mat, fmpz
    old = mp.mp.dps
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=delta, eta=0.51)
    mp.mp.dps = old
    return [[int(R[r, c]) for c in range(n + 1)] for r in range(n)]

def normalize(row):
    from math import gcd
    from functools import reduce
    coeffs = [int(c) for c in row[:-1]]
    g = reduce(gcd, [abs(c) for c in coeffs if c]) or 1
    coeffs = [c // g for c in coeffs]
    if next((c for c in coeffs if c), 1) < 0:
        coeffs = [-c for c in coeffs]
    return coeffs


# ---------------- corner quadrature with precision boost near g = 0 ----------------
def _mul(A, B, N):
    return [mp.fsum(A[i]*B[n - i] for i in range(n + 1)) for n in range(N + 1)]

def _sqrt(Q, N):
    b = [mp.sqrt(Q[0])]
    for n in range(1, N + 1):
        s = Q[n] - mp.fsum(b[k]*b[n - k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b

def integrand_series_hp(g, c, N):
    """u-series of the corner integrand at node g; internal precision boosted by
    proximity of g to 0 (fixes the (1+g)^(-1-u)*sqrt(...) - 1 cancellation)."""
    base = mp.mp.dps
    boost = 0
    if g < mp.mpf('0.1'):
        boost = int(-mp.log10(g)) + 15
    mp.mp.dps = base + boost
    g2 = +g
    c2 = +c
    L = mp.log(g2)
    m = mp.log1p(g2)
    fac = [mp.factorial(k) for k in range(N + 1)]
    E1 = [(-L)**k/fac[k] for k in range(N + 1)]
    E2 = [(-2*L)**k/fac[k] for k in range(N + 1)]
    Q = [1 + g2*g2 - 2*c2*g2] + [g2*g2*E2[k] - 2*c2*g2*E1[k] for k in range(1, N + 1)]
    B = _sqrt(Q, N)
    A = [(-m)**k/fac[k]/(1 + g2) for k in range(N + 1)]
    AB = _mul(A, B, N)
    AB[0] -= 1
    G = [L**k/fac[k]/g2 for k in range(N + 1)]
    out = _mul(AB, G, N)
    mp.mp.dps = base
    return out

def corner_coeffs_hp(a, N, pts=None):
    c = mp.cos(a)
    cache = {}
    def F(g, k):
        key = g
        if key not in cache:
            cache[key] = integrand_series_hp(g, c, N)
        return cache[key][k]
    if pts is None:
        pts = [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    return [mp.quad(lambda g: F(g, k), pts) for k in range(N + 1)]
