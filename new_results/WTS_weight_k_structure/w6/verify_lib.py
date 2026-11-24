"""verify_lib: independent quadrature + direct-definition check for the weight-6 (and nearby) claims.

Deliberately shares NO code with s02_eval.py (no mp.quad / tanh-sinh).  Differences:
  * composite Gauss-Legendre in y = -log(distance to the endpoint), nodes computed here by Newton on P_n;
  * log t and log(1-t) are evaluated from the exact factorisations
        t   = (v-b)(v+a)/(4v),     1-t = (a-v)(v+b)/(4v)
    as SUMS of logs, with v-b and a-v kept as exact distances (no log1p, no cancellation);
  * the reading of finals_K.json is re-implemented here.
PRECISION: mp.mp.dps is set explicitly (mp.mp.dps, never mp.dps) by every entry point AFTER all imports,
and this module imports nothing that resets it (no dilog_exact / m3_tails_exact / *_tails_exact.pkl).
"""
import os, json, math
from fractions import Fraction as Fr
import mpmath as mp

W = os.path.dirname(os.path.abspath(__file__))

POLE_LABEL = {((0, 0)): "0", ((1, 1)): "a", ((-1, -1)): "-a", ((-1, 1)): "b", ((1, -1)): "-b"}


def set_dps(n):
    mp.mp.dps = n          # mp.mp.dps, NOT mp.dps
    return mp.mp.dps


def fr2mp(f):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


def load_finals(path):
    """finals_K.json -> {(i,j): {label: (Fr p0, Fr p1)}} for the REAL simple-pole part, with
    assertions that there are no +-i poles, no constant, and every coefficient is real."""
    data = json.load(open(path))
    out = {}
    for ks, P in data["finals"].items():
        i, j = map(int, ks.split(","))
        assert not P["ipos"] and not P["ineg"], (ks, "i-poles")
        assert all(Fr(x) == 0 for x in P["const"]), (ks, "const")
        d = {}
        for r, c in P["real"]:
            assert Fr(r[2]) == 0 and Fr(r[3]) == 0 and Fr(c[2]) == 0 and Fr(c[3]) == 0
            lab = POLE_LABEL[(Fr(r[0]), Fr(r[1]))]
            assert lab not in d
            d[lab] = (Fr(c[0]), Fr(c[1]))
        out[(i, j)] = d
    return out


_GL_CACHE = {}


def gauss_legendre(n):
    """nodes/weights on [-1,1] at the current mp precision (Newton on P_n)."""
    key = (n, mp.mp.prec)
    if key in _GL_CACHE:
        return _GL_CACHE[key]
    xs, ws = [], []
    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for m in range(2, n + 1):
                p0, p1 = p1, ((2 * m - 1) * x * p1 - (m - 1) * p0) / m
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(2) ** (-mp.mp.prec + 8):
                break
        p0, p1 = mp.mpf(1), x
        for m in range(2, n + 1):
            p0, p1 = p1, ((2 * m - 1) * x * p1 - (m - 1) * p0) / m
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    _GL_CACHE[key] = (xs, ws)
    return xs, ws


def panel_edges(Y):
    """graded panel edges in y = -log(distance): the nearest singularities of the integrands in the y-plane are at
    y = -0.35 (log(2-u)) and y = -0.69 (t=1), so the first panels are short."""
    e = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0]
    while e[-1] < Y:
        e.append(e[-1] + 1.0)
    return e


def quad_pieces(finals, keys, dps, h=None, n=None, Y=None):
    """returns {key: value of int_b^a K_key(v) log^i t log^j(1-t) dv} with composite Gauss-Legendre.
    mp.mp.dps must already be `dps`.  Uses y = -log(distance to endpoint)."""
    assert mp.mp.dps == dps
    r2 = mp.sqrt(2)
    a = 1 + r2
    b = r2 - 1
    if n is None:
        n = 40 + int(dps * 0.36)          # 60 digits -> 61 ... 110 digits -> 79
    if Y is None:
        Y = math.log(10) * (dps + 14) + 30
    xs, ws = gauss_legendre(n)
    # pole data per key
    pdata = {}
    for key in keys:
        i, j = key
        lst = []
        for lab, (p0, p1) in finals[key].items():
            c = fr2mp(p0) + fr2mp(p1) * r2
            r = {"0": mp.mpf(0), "a": a, "-a": -a, "b": b, "-b": -b}[lab]
            lst.append((lab, r, c))
            if lab == "b":
                assert j >= 1, ("divergent", key)
            if lab == "a":
                assert i >= 1, ("divergent", key)
        pdata[key] = lst
    res = {key: mp.mpf(0) for key in keys}
    edges = panel_edges(Y)
    smaxL = 1 - b            # left half: v in (b, 1)
    smaxR = a - 1            # right half: v in (1, a)
    lsmaxL = mp.log(smaxL)
    lsmaxR = mp.log(smaxR)
    two = mp.mpf(2)
    for p in range(len(edges) - 1):
        y0, y1 = mp.mpf(edges[p]), mp.mpf(edges[p + 1])
        for xk, wk in zip(xs, ws):
            y = y0 + (y1 - y0) * (xk + 1) / 2
            wgt = wk * (y1 - y0) / 2
            e = mp.exp(-y)
            # ----- left half: v = b + s
            s = smaxL * e
            ls = lsmaxL - y
            v = b + s
            l4v = mp.log(4 * v)
            Lt = ls + mp.log(v + a) - l4v
            L1 = mp.log(two - s) + mp.log(v + b) - l4v
            jac = s * wgt
            for key in keys:
                i, j = key
                K = mp.mpf(0)
                for lab, r, c in pdata[key]:
                    K += c / (s if lab == "b" else (v - r))
                res[key] += jac * K * Lt ** i * L1 ** j
            # ----- right half: v = a - u
            u = smaxR * e
            lu = lsmaxR - y
            v = a - u
            l4v = mp.log(4 * v)
            Lt = mp.log(two - u) + mp.log(v + a) - l4v
            L1 = lu + mp.log(v + b) - l4v
            jac = u * wgt
            for key in keys:
                i, j = key
                K = mp.mpf(0)
                for lab, r, c in pdata[key]:
                    K += c / (-u if lab == "a" else (v - r))
                res[key] += jac * K * Lt ** i * L1 ** j
    return res


def weight_part(finals, d, dps, **kw):
    """W_{d+1} = sum over keys i+j=d, i>=j of the piece integrals."""
    keys = [(i, d - i) for i in range(d, (d - 1) // 2, -1) if (i, d - i) in finals]
    pcs = quad_pieces(finals, keys, dps, **kw)
    return sum(pcs.values(), mp.mpf(0)), pcs


# ------------------------------------------------------------------------------------------
# direct definition: M_k = -F^(k)(0)/(2 k!) = -(1/2) int_0^1 [eps^k] Phi dt,
# Phi = A e^{eps l1} + B e^{eps l2} - sqrt(A^2 e^{2 eps l1} + B^2 e^{2 eps l2}),  A=1/t, B=1/(1-t).
# Here sqrt(...) = Q' sqrt(1+u) with Q' = q/(t(1-t)), q^2 = t^2+(1-t)^2, u = alpha X + beta Y.
# Phi(t) is symmetric under t <-> 1-t, so M_k = -int_0^{1/2} [eps^k] Phi dt.
# ------------------------------------------------------------------------------------------

def _series_mul(P, Q, k):
    out = [mp.mpf(0)] * (k + 1)
    for i, pi_ in enumerate(P):
        if pi_ == 0:
            continue
        for j, qj in enumerate(Q):
            if i + j > k:
                break
            out[i + j] += pi_ * qj
    return out


def phi_coeff(t, one_minus_t, l1, l2, k):
    """[eps^k] Phi as a function of t, computed from truncated eps-series arithmetic.
    one_minus_t passed separately to avoid cancellation errors for t near 1 (we only use t<=1/2 anyway)."""
    t2 = t * t
    s2 = one_minus_t * one_minus_t
    q2 = t2 + s2
    q = mp.sqrt(q2)
    alpha = s2 / q2
    beta = t2 / q2
    fact = [mp.factorial(m) for m in range(k + 1)]
    X = [mp.mpf(0)] + [(2 * l1) ** m / fact[m] for m in range(1, k + 1)]
    Yv = [mp.mpf(0)] + [(2 * l2) ** m / fact[m] for m in range(1, k + 1)]
    u = [alpha * x + beta * y for x, y in zip(X, Yv)]
    # sqrt(1+u) = sum_j binom(1/2, j) u^j, u has zero constant term
    tot = [mp.mpf(0)] * (k + 1)
    tot[0] = mp.mpf(1)
    pw = [mp.mpf(1)] + [mp.mpf(0)] * k
    bj = mp.mpf(1)
    for j in range(1, k + 1):
        bj = bj * (mp.mpf(1) / 2 - (j - 1)) / j
        pw = _series_mul(pw, u, k)
        for m in range(k + 1):
            tot[m] += bj * pw[m]
    Qp = q / (t * one_minus_t)
    return l1 ** k / fact[k] / t + l2 ** k / fact[k] / one_minus_t - Qp * tot[k]


def direct_M(k, dps, h=None, n=None, Y=None):
    """M_k from the definition.  mp.mp.dps must be the WORKING precision (>= target + ~0.45*Y);
    returns M_k truncated at t_min = exp(-Y)/2.  Returns (value, Y)."""
    if n is None:
        n = 48 + int(dps * 0.15)
    xs, ws = gauss_legendre(n)
    edges = panel_edges(Y)
    half = mp.mpf(1) / 2
    tot = mp.mpf(0)
    for p in range(len(edges) - 1):
        y0, y1 = mp.mpf(edges[p]), mp.mpf(edges[p + 1])
        for xk, wk in zip(xs, ws):
            y = y0 + (y1 - y0) * (xk + 1) / 2
            t = half * mp.exp(-y)
            l1 = -mp.log(2) - y
            omt = 1 - t
            l2 = mp.log(omt)
            tot += wk * (y1 - y0) / 2 * t * phi_coeff(t, omt, l1, l2, k)
    return -tot
