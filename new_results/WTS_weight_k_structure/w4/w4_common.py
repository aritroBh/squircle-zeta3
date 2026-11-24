"""Shared helpers for the weight-4 study.  PRECISION: set mp.mp.dps AFTER all imports (call set_dps)."""
import sys, os, pickle, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import mpmath as mp
import fastred                       # exact Q(sqrt2,i) arithmetic; imports sympy only (does not touch mp.mp.dps)
from fractions import Fraction as Fr

RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)

def set_dps(n):
    mp.mp.dps = n                    # NOT mp.dps = n (that does nothing)
    return mp.mp.dps

def frac_mp(f):
    f = Fr(f); return mp.mpf(f.numerator) / mp.mpf(f.denominator)

def qsi_real_mp(q):
    """QSI (p0 + p1 sqrt2 + i(p2 + p3 sqrt2)) -> mpmath real; asserts the imaginary part is exactly 0."""
    assert q.p2 == 0 and q.p3 == 0, ("imaginary coefficient", q)
    return frac_mp(q.p0) + frac_mp(q.p1) * mp.sqrt(2)

def load_red(k):
    return pickle.load(open(os.path.join(HERE, f"red_{k}.pkl"), "rb"))

def pieces_mp(k, d):
    """final pieces of weight d+1 (i+j = d): list of (i, j, [(pole_mp, coef_mp), ...]).
    Asserts: simple real poles only, no i-poles, no constant."""
    red = load_red(k)
    out = []
    for (i, j), P in sorted(red.items(), key=lambda kv: (-kv[0][0], kv[0][1])):
        if i + j != d: continue
        real, ipos, ineg, const = P
        assert not ipos and not ineg and const.is_zero(), ("non-final piece", k, (i, j))
        out.append((i, j, [(qsi_real_mp(r), qsi_real_mp(c)) for r, c in real.items()]))
    return out

def make_integrand(pieces):
    """f(v) = sum K_ij(v) log^i t log^j(1-t) with t = (v-b)(v+a)/(4v), 1-t = (a-v)(v+b)/(4v).
    Logs are formed from the factored form so no cancellation occurs next to v = b (t->0) or v = a (t->1)."""
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
    def f(v):
        l4v = mp.log(4 * v)
        L1 = mp.log(v - b) + mp.log(v + a) - l4v
        L2 = mp.log(a - v) + mp.log(v + b) - l4v
        tot = mp.mpf(0)
        for (i, j, poles) in pieces:
            K = mp.mpf(0)
            for r, c in poles:
                K += c / (v - r)
            tot += K * L1**i * L2**j
        return tot
    return f

def integrate(pieces, label=""):
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
    f = make_integrand(pieces)
    val, err = mp.quad(f, [b, 1, a], error=True)
    return val, err
