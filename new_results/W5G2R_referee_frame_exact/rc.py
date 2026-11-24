"""W5G2R referee: own primitives. Own word enumeration, own Galois action, own shuffle, own numerical evaluator.
Project files are read only as DATA (pickled rows). No project module is imported here except that unpickling
validated_w5.pkl needs the qsqrt2 class (provenance field only)."""
import sys, os, pickle, time, itertools
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
ROOT = '<project>'
T = ROOT + '/research/08_theory'
HERE = T + '/W5G2R_referee_frame_exact'
sys.path.append(ROOT + '/paper/squircle_note/scripts')   # only for unpickling Q(sqrt2) provenance objects

def log(s, f):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    with open(os.path.join(HERE, f), 'a') as fh: fh.write(s + '\n')

# letters: -1 = dt/t ; j in 0..7 = z^j dt/(1 - z^j t), z = exp(2 pi i/8). H(a1..an) = int_{1>t1>..>tn>0}.
# admissible: a1 != 0, an != -1
@lru_cache(None)
def WORDS(n):
    if n == 0: return ((),)
    return tuple(w for w in itertools.product(range(-1, 8), repeat=n) if w[0] != 0 and w[-1] != -1)

def act(w, u): return tuple(x if x == -1 else (u * x) % 8 for x in w)
def actv(d, u): return {act(w, u): c for w, c in d.items()}

@lru_cache(None)
def shw(a, b):
    """shuffle of two words -> dict"""
    if not a: return {b: 1}
    if not b: return {a: 1}
    out = {}
    for w, c in shw(a[1:], b).items():
        k = (a[0],) + w; out[k] = out.get(k, 0) + c
    for w, c in shw(a, b[1:]).items():
        k = (b[0],) + w; out[k] = out.get(k, 0) + c
    return out

def shv(x, y):
    out = {}
    for a, ca in x.items():
        for b, cb in y.items():
            for w, c in shw(a, b).items(): out[w] = out.get(w, 0) + ca * cb * c
    return {w: c for w, c in out.items() if c}

def prim(d):
    """primitive integer form (sign: first nonzero in sorted word order positive)"""
    from math import gcd, lcm
    d = {w: F(c) for w, c in d.items() if c}
    if not d: return {}
    L = lcm(*[c.denominator for c in d.values()]); items = sorted(d.items())
    ints = [(w, int(c * L)) for w, c in items]; g = 0
    for _, a in ints: g = gcd(g, abs(a))
    if ints[0][1] < 0: g = -g
    return {w: a // g for w, a in ints}

E = {(2,): 1, (6,): -1}       # e = H(2) - H(6), value i pi/2

def load_rows(fam):
    """W5C1 literal rows (data) -> list of (label, dict word->int) with words decoded by my own enumeration"""
    W5 = WORDS(5)
    R = pickle.load(open(T + f'/W5C1_certificate_scoping/rels_{fam}.pkl', 'rb'))
    return [(lab, {W5[i]: c for i, c in row}) for lab, row in R]

def bases():
    out = {}
    for n in (1, 2, 3):
        sel = pickle.load(open(T + f'/A3c_claude_elementary/Q/weight{n}.pickle', 'rb'))[2]
        out[n] = [{tuple(w): int(c) for w, c in row} for row in sel]
    T0 = pickle.load(open(T + '/A3h_weight4_target/w4_systems_T.pkl', 'rb'))['T0']
    S = pickle.load(open(T + '/A3h5_exact_certificate/stage2_rows_S.pkl', 'rb'))['S']
    b4 = []
    for i in S:
        d = {}
        for w, (q, c) in T0[i][1]:
            assert q == 0; d[tuple(w)] = F(c)
        b4.append(d)
    out[4] = b4
    return out

# ---------------- numerical evaluator ----------------
def _gl(q):
    x, w = np.polynomial.legendre.leggauss(q)
    # integration matrix on [-1,1]: S[i,j] = int_{-1}^{x_i} l_j(s) ds, via Legendre expansion
    V = np.polynomial.legendre.legvander(x, q - 1)            # V[i,k] = P_k(x_i)
    Vi = np.linalg.inv(V)                                       # coefficients of l_j: c[k,j] = Vi[k,j]
    # antiderivative of P_k from -1: (P_{k+1} - P_{k-1})/(2k+1), k>=1; P_0 -> x+1
    A = np.zeros((q, q))
    for k in range(q):
        c = np.zeros(q + 1); c[k] = 1
        ci = np.polynomial.legendre.legint(c, lbnd=-1)
        A[:, k] = np.polynomial.legendre.legval(x, ci)
    return x, w, A @ Vi

def evaluate(nmax=5, q=24, K=62, logf=None):
    """values H(w) for every word w of length 1..nmax with last letter != -1 (first letter arbitrary),
    returned as dict for admissible words (first != 0). Cumulative Gauss-Legendre on dyadic panels toward 0 and 1."""
    z = np.exp(2j * np.pi / 8)
    x, wq, Sm = _gl(q)
    # suffix-closed list of words by length (last letter != -1)
    words = {0: [()]}
    for L in range(1, nmax + 1):
        words[L] = [w for w in itertools.product(range(-1, 8), repeat=L) if w[-1] != -1]
    idx = {L: {w: i for i, w in enumerate(words[L])} for L in words}
    # child index: word a.w' -> (a, index of w')
    child = {L: (np.array([w[0] for w in words[L]]), np.array([idx[L - 1][w[1:]] for w in words[L]])) for L in range(1, nmax + 1)}
    Fv = {L: np.zeros(len(words[L]), complex) for L in words}; Fv[0][:] = 1
    panels = [('t', 0.0, 2.0 ** -K)] + [('t', 2.0 ** -(k + 1), 2.0 ** -k) for k in range(K - 1, 0, -1)] \
           + [('s', 2.0 ** -k, 2.0 ** -(k + 1)) for k in range(1, K)] + [('s', 2.0 ** -K, 0.0)]
    for kind, a, b in panels:
        h = abs(b - a) / 2
        if kind == 't':
            t = (a + b) / 2 + h * x; s = 1 - t
        else:      # a > b in s; t increases as s decreases
            s = (a + b) / 2 - h * x; t = 1 - s
        # forms omega_a(t) for a = -1..7 at nodes
        om = np.zeros((9, q), complex)
        om[0] = 1 / t
        om[1] = 1 / s
        for j in range(1, 8):
            zj = z ** j
            om[j + 1] = zj / ((1 - zj) + zj * s)
        Fn = {0: np.ones((1, q), complex)}
        newF = {0: Fv[0]}
        for L in range(1, nmax + 1):
            let, ch = child[L]
            g = om[let + 1] * Fn[L - 1][ch]                     # (nw, q)
            Fn[L] = Fv[L][:, None] + h * (g @ Sm.T)
            newF[L] = Fv[L] + h * (g @ wq)
        Fv = newF
    out = {}
    for L in range(1, nmax + 1):
        for w, v in zip(words[L], Fv[L]):
            if w[0] != 0: out[w] = v
    return out

def val(d, H):
    """value of a row and its scale sum |c||H|"""
    s = 0j; sc = 0.0
    for w, c in d.items():
        h = H[w]; s += float(c) * h if abs(c) < 2**52 else complex(F(c) * F(h.real)) ; sc += abs(float(c)) * abs(h)
    return s, sc
