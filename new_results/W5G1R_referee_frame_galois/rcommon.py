"""Referee helpers (W5G1R). Independent of W5G1 code: own Galois action, own shuffle, own numerical evaluator.
Conventions (read from A3c_claude_elementary/Q/cmzv.py): letter -1 = dt/t, letter j in 0..7 = z^j dt/(1 - z^j t),
z = exp(2 pi i/8); word (a1..an) = iterated integral over 1 > t1 > ... > tn > 0 of prod omega_{ak}(tk).
Admissible: a1 != 0, an != -1."""
import sys, time, os
sys.dont_write_bytecode = True
from fractions import Fraction as F
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = '<project>'
T = ROOT + '/research/08_theory'
def log(s, f):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    with open(os.path.join(HERE, f), 'a') as fh: fh.write(s + '\n')
def act(w, u): return tuple(x if x == -1 else (u * x) % 8 for x in w)
def actv(d, u): return {act(w, u): c for w, c in d.items()}
_shc = {}
def sh(a, b):
    k = (a, b)
    if k in _shc: return _shc[k]
    if not a: r = {b: 1}
    elif not b: r = {a: 1}
    else:
        r = {}
        for w, c in sh(a[1:], b).items(): r[(a[0],) + w] = r.get((a[0],) + w, 0) + c
        for w, c in sh(a, b[1:]).items(): r[(b[0],) + w] = r.get((b[0],) + w, 0) + c
    _shc[k] = r
    return r
def shv(x, y):
    out = {}
    for a, ca in x.items():
        for b, cb in y.items():
            for w, c in sh(a, b).items(): out[w] = out.get(w, 0) + ca * cb * c
    return {w: c for w, c in out.items() if c}
def words(n):
    from itertools import product
    if n == 0: return [()]
    return [w for w in product(range(-1, 8), repeat=n) if w[0] != 0 and w[-1] != -1]

# ---------- independent double-precision evaluator: graded panels + Chebyshev cumulative integration ----------
_nodes = None
def _grid(K=48, m=24):
    global _nodes
    if _nodes is not None: return _nodes
    bps = [0.0] + [2.0 ** -k for k in range(K, 1, -1)] + [0.5] + [1 - 2.0 ** -k for k in range(2, K + 1)]
    x = np.cos(np.pi * (np.arange(m) + 0.5) / m)[::-1]   # Chebyshev-Gauss nodes on (-1,1)
    # cumulative integration matrix on (-1,1): values at nodes -> integral from -1 to node
    V = np.polynomial.chebyshev.chebvander(x, m - 1)
    Vi = np.linalg.inv(V)
    Ci = np.zeros((m, m))
    for k in range(m):
        e = np.zeros(m); e[k] = 1
        ci = np.polynomial.chebyshev.chebint(e, lbnd=-1)
        Ci[:, k] = np.polynomial.chebyshev.chebval(x, ci)
    Mcum = Ci @ Vi                                   # cumulative from -1
    wtot = np.array([np.polynomial.chebyshev.chebval(1, np.polynomial.chebyshev.chebint(np.eye(m)[k], lbnd=-1)) for k in range(m)]) @ Vi
    pans = []
    for a, b in zip(bps[:-1], bps[1:]):
        pans.append((a, b, a + (b - a) * (x + 1) / 2))
    _nodes = (pans, Mcum, wtot)
    return _nodes
Z8 = np.exp(2j * np.pi * np.arange(8) / 8)
Z8[0], Z8[2], Z8[4], Z8[6] = 1, 1j, -1, -1j
def _form(a, t):
    if a == -1: return 1.0 / t
    return Z8[a] / (1 - Z8[a] * t)
_hc = {}
def H(w):
    """H(w) at t=1 for an admissible word, double precision (~1e-12)."""
    w = tuple(w)
    if w in _hc: return _hc[w]
    pans, Mcum, wtot = _grid()
    # g_k(t) = iterated integral of word w[k:] from 0 to t ; build from the inside out, panel by panel
    g = [np.ones(len(p[2]), dtype=complex) for p in pans]      # word ()
    for a in reversed(w[1:]):
        acc = 0j; new = []
        for (lo, hi, t), gv in zip(pans, g):
            f = _form(a, t) * gv
            cum = acc + (hi - lo) / 2 * (Mcum @ f)
            new.append(cum)
            acc = acc + (hi - lo) / 2 * (wtot @ f)
        g = new
    a = w[0]; tot = 0j
    for (lo, hi, t), gv in zip(pans, g):
        tot += (hi - lo) / 2 * (wtot @ (_form(a, t) * gv))
    _hc[w] = tot
    return tot
def val(d): return sum(complex(float(c), 0) * H(w) for w, c in d.items() if c)
