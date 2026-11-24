"""verify_lib.py -- my OWN quadrature and helpers for the adversarial check of the weight-4 results.
Nothing here imports w4_common / w4_eval / w4_known / id_lib; the only thing taken from the author's folder is
fastred.py (needed solely so that the exact pickles red_<k>.pkl can be unpickled) and the pickles themselves.

Quadrature: Gauss-Legendre panels in x = -log(t) (lower half v in [b,1], t in (0,1/2]) and in y = -log(1-t)
(upper half v in [1,a]).  In these variables log t / log(1-t) are exact (one of them is -x), the integrand is analytic on the
real axis, the end-point log singularities are gone, and v-b, a-v are formed from the exact factorised forms
  v - b = 4 v t/(v+a),   a - v = 4 v (1-t)/(v+b)
so there is no cancellation next to the end points.   PRECISION: call set_dps(n) AFTER all imports; it sets mp.mp.dps.
"""
import sys, os, pickle, json, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp

W4 = "<project>/research/08_theory/WTS_weight_k_structure/w4"
if W4 not in sys.path:
    sys.path.insert(0, W4)
import fastred                              # exact QSI arithmetic (needed to unpickle only)

POLE_NAMES = {(0, 0): "0", (1, 1): "a", (-1, -1): "ma", (-1, 1): "b", (1, -1): "mb"}   # (p0,p1) of p0+p1*sqrt2


def set_dps(n):
    mp.mp.dps = n                           # NOT mp.dps = n
    return mp.mp.dps


def fr_mp(f):
    f = Fr(f)
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


def qsi_mpc(q):
    r2 = mp.sqrt(2)
    return mp.mpc(fr_mp(q.p0) + fr_mp(q.p1) * r2, fr_mp(q.p2) + fr_mp(q.p3) * r2)


def red_exact(k, folder=None):
    """exact pieces; red_<k>.pkl from the author's folder, or (k>10) from the directory in $VERIFY_EXTRA (my own re-run)"""
    if folder is None:
        folder = W4
        if not os.path.exists(os.path.join(W4, f"red_{k}.pkl")):
            folder = os.environ["VERIFY_EXTRA"]
    return pickle.load(open(os.path.join(folder, f"red_{k}.pkl"), "rb"))


def pole_name(r):
    assert r.p2 == 0 and r.p3 == 0, ("non-real pole", r)
    return POLE_NAMES[(r.p0, r.p1)]


def pieces_mp(red):
    """exact red dict -> {(i,j): dict(real={pole: mpc}, ipos=[..], ineg=[..], const=mpc)} at the current precision"""
    out = {}
    for key, P in red.items():
        real, ipos, ineg, const = P
        out[key] = dict(real={pole_name(r): qsi_mpc(c) for r, c in real.items()},
                        ipos=[qsi_mpc(e) for e in ipos], ineg=[qsi_mpc(f) for f in ineg], const=qsi_mpc(const))
    return out


# ---------------------------------------------------------------- Gauss-Legendre
_RULES = {}


def gl_rule(n):
    key = (n, mp.mp.prec)
    if key in _RULES:
        return _RULES[key]
    xs = []
    ws = []
    half = (n + 1) // 2
    tol = mp.mpf(2) ** (-(mp.mp.prec - 10))
    for i in range(1, half + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(200):
            p0, p1 = mp.mpf(1), x
            for kk in range(2, n + 1):
                p0, p1 = p1, ((2 * kk - 1) * x * p1 - (kk - 1) * p0) / kk
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < tol:
                break
        p0, p1 = mp.mpf(1), x
        for kk in range(2, n + 1):
            p0, p1 = p1, ((2 * kk - 1) * x * p1 - (kk - 1) * p0) / kk
        dp = n * (x * p1 - p0) / (x * x - 1)
        w = 2 / ((1 - x * x) * dp * dp)
        xs.append(x)
        ws.append(w)
    nodes = []
    for x, w in zip(xs, ws):
        nodes.append((x, w))
        if not (n % 2 == 1 and x == xs[-1]):
            nodes.append((-x, w))
    if n % 2 == 1:                          # middle node was computed once; make sure it is exactly 0 weight-wise ok
        pass
    assert len(nodes) == n, (len(nodes), n)
    _RULES[key] = nodes
    return nodes


# ---------------------------------------------------------------- nodes
_CONST = {}


def consts():
    key = mp.mp.prec
    if key not in _CONST:
        r2 = mp.sqrt(2)
        _CONST[key] = (r2, 1 + r2, r2 - 1, mp.log(2))
    return _CONST[key]


def make_node(x, half):
    """half 0: lower half v in [b,1], t = exp(-x).  half 1: upper half v in [1,a], 1-t = exp(-x)."""
    r2, a, b, _ = consts()
    e = mp.exp(-x)
    if half == 0:
        t = e
        v = (2 * t - 1) + mp.sqrt(2 - 4 * t + 4 * t * t)
        l1 = -x
        l2 = mp.log1p(-t)
        dif = {"0": v, "a": v - a, "ma": v + a, "b": 4 * v * t / (v + a), "mb": v + b}
    else:
        s = e
        w = (2 * s - 1) + mp.sqrt(2 - 4 * s + 4 * s * s)
        v = 1 / w
        l2 = -x
        l1 = mp.log1p(-s)
        dif = {"0": v, "a": -4 * v * s / (v + b), "ma": v + a, "b": v - b, "mb": v + b}
    meas = 4 * v * v / (v * v + 1) * e      # dv/dt * |dt/dx|
    return dict(x=x, half=half, v=v, l1=l1, l2=l2, dif=dif, meas=meas, e=e)


def eval_pf(P, nd):
    """value of the partial-fraction function P at the node (real simple poles via the stable differences)."""
    v = nd["v"]
    tot = mp.mpc(0) if (P["ipos"] or P["ineg"] or P["const"] != 0 or any(c.imag != 0 for c in P["real"].values())) else mp.mpf(0)
    for name, c in P["real"].items():
        tot += (c if isinstance(tot, mp.mpc) else c.real) / nd["dif"][name]
    if P["ipos"] or P["ineg"] or P["const"] != 0:
        for j, c in enumerate(P["ipos"], start=1):
            tot += c / mp.mpc(v, -1) ** j
        for j, c in enumerate(P["ineg"], start=1):
            tot += c / mp.mpc(v, 1) ** j
        tot += P["const"]
    return tot


def piece_integrand(pieces, keys):
    """callable nd -> sum_{key in keys} K_key(v) l1^i l2^j * dv/dt*|dt/dx|"""
    def f(nd):
        tot = 0
        for key in keys:
            i, j = key
            tot += eval_pf(pieces[key], nd) * nd["l1"] ** i * nd["l2"] ** j
        return tot * nd["meas"]
    return f


def panel_list(X, wmax=6):
    x0 = mp.log(2)
    widths = [mp.mpf(3) / 10, mp.mpf(3) / 10, mp.mpf(2) / 5, mp.mpf(1) / 2, mp.mpf(3) / 4, 1, mp.mpf(3) / 2, 2, 3, 4, 5]
    edges = [x0]
    idx = 0
    while edges[-1] < X:
        w = widths[idx] if idx < len(widths) else wmax
        edges.append(edges[-1] + w)
        idx += 1
    return list(zip(edges[:-1], edges[1:]))


def integrate(fs, X, n, halves=(0, 1), wmax=6):
    """fs: list of callables f(nd) (full integrand per dx incl. measure).  Returns list of integrals over x in [log2, X]
    (summed over the requested halves)."""
    rule = gl_rule(n)
    tot = [mp.mpf(0) for _ in fs]
    for half in halves:
        for lo, hi in panel_list(X, wmax):
            c = (lo + hi) / 2
            h = (hi - lo) / 2
            for xi, wi in rule:
                nd = make_node(c + h * xi, half)
                for idx, f in enumerate(fs):
                    tot[idx] += h * wi * f(nd)
    return tot


def default_X(D):
    return (D + 12) * mp.log(10) + 12


# ---------------------------------------------------------------- direct M_k from the squircle definition
def mk_direct_integrand(k):
    """G(x) e^{-x}, G(x) = [eps^k] Phi at t = e^{-x} in (0,1/2], Phi = t^{eps-1} + (1-t)^{eps-1} - sqrt(t^{2eps-2}+(1-t)^{2eps-2}).
    M_k = -(1/2) int_0^1 [eps^k]Phi dt = - int_0^{1/2} [eps^k]Phi dt  (Phi is symmetric under t -> 1-t)."""
    fk = mp.factorial
    def f(nd):
        assert nd["half"] == 0
        x = nd["x"]
        t = nd["e"]
        omt = -mp.expm1(-x)
        l1 = nd["l1"]
        l2 = nd["l2"]
        A = 1 / t
        B = 1 / omt
        r0 = (t / omt) ** 2
        delta = l2 - l1
        # r(eps) = r0 exp(2 eps delta)
        r = [r0 * (2 * delta) ** n / fk(n) for n in range(k + 1)]
        s = [mp.sqrt(1 + r0)]
        for n in range(1, k + 1):
            acc = r[n] - sum(s[m] * s[n - m] for m in range(1, n))
            s.append(acc / (2 * s[0]))
        u = [-r0 / (1 + s[0])] + [-s[n] for n in range(1, k + 1)]      # (1 - sqrt(1+r))_n
        coef = A * sum(u[n] * l1 ** (k - n) / fk(k - n) for n in range(k + 1)) + B * l2 ** k / fk(k)
        return -coef * nd["e"]                                           # M_k = - int coef dt, dt = e^{-x} dx
    return f


# ---------------------------------------------------------------- exact k-table of the weight-4 constants
def exact_table():
    """{(key, pole): [ (p0,p1) for k^0, k^1, k^2 ]} with coefficient of k^j = p0 + p1 sqrt2 (Fractions),
    from the exact quadratic through k=4,5,6 (the caller tests it at other k)."""
    T = {}
    for kk in (4, 5, 6):
        red = red_exact(kk)
        for key in ((3, 0), (2, 1)):
            real, ip, im, const = red[key]
            assert not ip and not im and const.is_zero()
            for r, c in real.items():
                T.setdefault((key, pole_name(r)), {})[kk] = (Fr(c.p0), Fr(c.p1))
    out = {}
    for kp, d in T.items():
        co = []
        for comp in (0, 1):
            ys = [d.get(kk, (Fr(0), Fr(0)))[comp] for kk in (4, 5, 6)]
            c2 = ((ys[2] - ys[1]) - (ys[1] - ys[0])) / 2
            c1 = (ys[1] - ys[0]) - c2 * 9
            c0 = ys[0] - c1 * 4 - c2 * 16
            co.append((c0, c1, c2))
        out[kp] = [(co[0][j], co[1][j]) for j in range(3)]
    return out


def k3_table():
    red = red_exact(3)
    out = {}
    for key in ((3, 0), (2, 1)):
        real, ip, im, const = red[key]
        assert not ip and not im and const.is_zero()
        for r, c in real.items():
            out[(key, pole_name(r))] = (Fr(c.p0), Fr(c.p1))
    return out


def piece_dict(coefs):
    """coefs {(key,pole): (p0,p1)} Fractions -> pieces dict usable by piece_integrand (current precision)"""
    r2 = mp.sqrt(2)
    out = {}
    for (key, pole), (p0, p1) in coefs.items():
        d = out.setdefault(key, dict(real={}, ipos=[], ineg=[], const=mp.mpc(0)))
        d["real"][pole] = mp.mpc(fr_mp(p0) + fr_mp(p1) * r2)
    return out
