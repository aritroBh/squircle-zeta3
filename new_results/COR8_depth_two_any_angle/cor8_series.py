"""cor8_series.py -- Puiseux series in c and log c with mpf coefficients.
A series is dict (k, l) -> mpf, meaning sum c^k (log c)^l, k >= -K0, l >= 0.
Working precision set by caller. Truncation: k <= K."""
import mpmath as mp

K = 60          # max c-power kept
MAXL = 4

def trim(S):
    return {(k, l): v for (k, l), v in S.items() if v and k <= K and l <= MAXL}
def sadd(*args):
    out = {}
    for S in args:
        for k, v in S.items():
            out[k] = out.get(k, 0) + v
    return trim(out)
def sscale(S, a):
    return trim({k: v*a for k, v in S.items()})
def smul(A, B):
    out = {}
    for (k1, l1), v1 in A.items():
        for (k2, l2), v2 in B.items():
            key = (k1+k2, l1+l2)
            out[key] = out.get(key, 0) + v1*v2
    return trim(out)
def sshift(S, dk):
    return trim({(k+dk, l): v for (k, l), v in S.items()})
def sintegrate(S):
    """termwise integral of c^k (log c)^l dc (k != -1)."""
    out = {}
    for (k, l), v in S.items():
        assert k != -1
        # int c^k log^l = c^{k+1} sum_{j=0}^l (-1)^j l!/(l-j)! log^{l-j} / (k+1)^{j+1}
        for j in range(0, l+1):
            import math
            coef = ((-1)**j) * math.factorial(l)//math.factorial(l-j) / mp.mpf(k+1)**(j+1)
            out[(k+1, l-j)] = out.get((k+1, l-j), 0) + v*coef
    return trim(out)
def sdiff(S):
    out = {}
    for (k, l), v in S.items():
        if k: out[(k-1, l)] = out.get((k-1, l), 0) + v*k
        if l: out[(k-1, l-1)] = out.get((k-1, l-1), 0) + v*l
    return trim(out)
def sexp_pow(base, n):
    """series^n (integer n >= 0), base has only nonnegative k? (not required if no (0,0) issues)"""
    out = {(0, 0): mp.mpf(1)}
    for _ in range(n):
        out = smul(out, base)
    return out
def srat_series(expr, c):
    """Taylor/Puiseux series of a rational function expr in symbol c (mpf coeffs)."""
    import sympy as sp
    num, den = sp.fraction(sp.together(expr))
    num = sp.Poly(sp.expand(num), c)
    den = sp.Poly(sp.expand(den), c)
    nd = {i: mp.mpf(sp.Rational(co)) for i, co in enumerate(num.all_coeffs()[::-1])}
    dd = {i: mp.mpf(sp.Rational(co)) for i, co in enumerate(den.all_coeffs()[::-1])}
    dmin = min(dd)
    # normalize denominator to (1 + u) with u having positive powers
    d0 = dd[dmin]
    u = {(i - dmin, 0): co/d0 for i, co in dd.items() if i != dmin}
    # 1/(1+u) = sum (-u)^m
    inv = {(0, 0): mp.mpf(1)}
    um = {(0, 0): mp.mpf(1)}
    uu = dict(u)
    for m in range(1, K+2):
        um = smul(um, uu)
        inv = sadd(inv, sscale(um, (-1)**m))
    inv = sscale(inv, 1/d0)
    out = {(0, 0): mp.mpf(0)}
    for i, co in nd.items():
        out = sadd(out, sscale(sshift(inv, i), co))
    return sshift(out, -dmin)
def _log1p(u):
    """log(1+u) for series u with positive k-powers only."""
    out = {}
    um = dict(u)
    for m in range(1, K+2):
        out = sadd(out, sscale(um, (-1)**(m+1)/mp.mpf(m)))
        um = smul(um, u)
    return trim(out)
def slog_of_rat(expr, c):
    """series of log(rat): rat = (nc[nmin]/dc[dmin]) * c^(nmin-dmin) * (1+u_num)/(1+u_den),
    = log|g0| + (nmin-dmin) log c + log(1+u_num) - log(1+u_den)."""
    import sympy as sp
    num, den = sp.fraction(sp.together(expr))
    num = sp.Poly(sp.expand(num), c); den = sp.Poly(sp.expand(den), c)
    nc = num.all_coeffs()[::-1]; dc = den.all_coeffs()[::-1]
    nmin = next(i for i, co in enumerate(nc) if co); dmin = next(i for i, co in enumerate(dc) if co)
    g0 = abs(mp.mpf(sp.Rational(nc[nmin])) / mp.mpf(sp.Rational(dc[dmin])))
    out = {(0, 0): mp.log(g0)}
    m = nmin - dmin
    if m:
        out[(0, 1)] = out.get((0, 1), 0) + m
    u_num = {(i - nmin, 0): mp.mpf(sp.Rational(co))/mp.mpf(sp.Rational(nc[nmin]))
             for i, co in enumerate(nc) if i > nmin}
    u_den = {(i - dmin, 0): mp.mpf(sp.Rational(co))/mp.mpf(sp.Rational(dc[dmin]))
             for i, co in enumerate(dc) if i > dmin}
    if u_num:
        out = sadd(out, _log1p(u_num))
    if u_den:
        out = sadd(out, sscale(_log1p(u_den), -1))
    return trim(out)
def sLi(n, expr, c):
    """series of Li_n(g(c)) for rational g with g(0) != 0,1 finite; or g(0)=0.
    Uses dLi_n/dc = Li_{n-1} g'/g, constants Li_n(g0) (mpf)."""
    import sympy as sp
    g = sp.sympify(expr) if isinstance(expr, str) else expr
    g0 = g.subs(c, 0)
    if g0 == 0:
        const = {1: mp.mpf(0)}
        cur = sscale(slog_of_rat(1 - g, c), -1)   # Li_1 = -log(1-g)
    else:
        g0m = mp.mpf(g0)
        const = {1: mp.re(mp.polylog(1, g0m))}
        cur = {(0, 0): const[1]}
    logg = slog_of_rat(g, c)
    dlogg = sdiff(logg)
    for k in range(2, n+1):
        cur = sintegrate(smul(cur, dlogg))
        if g0 == 0:
            const[k] = mp.mpf(0)
        else:
            const[k] = mp.re(mp.polylog(k, mp.mpf(g0)))
        cur[(0, 0)] = const[k]
    return trim(cur)
def sLi_inv(n, expr, c):
    """series of Li_n(g(c)) where g has a pole at 0: use inversion Li_n(g) = -(-1)^n Li_n(1/g) + log-poly.
    For n=2: Li2(g)+Li2(1/g) = -pi^2/6 - 1/2 log(-g)^2  (g not in (0,1]).
    For n=3: Li3(g)-Li3(1/g) = -pi^2/6 log(-g) - 1/6 log(-g)^3.
    For n=4: Li4(g)+Li4(1/g) = -pi^4/45 - pi^2/12 log(-g)^2 - 1/24 log(-g)^4.
    Real parts with log(-g) = log|g| + i*pi*sgn: we take the real unipotent convention;
    for real g>1 or g<0 with the i0-prescription, the real value uses log|g| and the
    formulas below verbatim with log(-g) -> log|g| for the REAL parts (even powers)
    -- handled by caller via slog_of_rat(|g|)."""
    raise NotImplementedError
