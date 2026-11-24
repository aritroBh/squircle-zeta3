"""W4RX core: an independent implementation (imports no module of the notes).

Words.  A level-8 word is a tuple of letters; letter m in 0..7 is the point zeta^m,
zeta = exp(2 pi i/8); letter 8 (ZERO) is the point 0.  The word stands for
    G(a1,...,an;1),  G(a1,...,an;y) = int_0^y dt/(t-a1) G(a2,...,an;t)   (a1 outermost).
Admissible = convergent: a1 != 1 (letter 0) and an != 0 (letter 8).
A vector is a dict word -> coefficient (Fraction, or K for Q(sqrt2) coefficients).

The notes use H-words (letters -1 = dt/t, j = zeta^j dt/(1-zeta^j t)); the dictionary is
    H(-1)-letter <-> G-letter 0 (8),   H-letter j <-> G-letter zeta^(-j) with a factor -1,
so H(w) = (-1)^{#letters != -1} G(g(w)).  It is used ONLY to compare with the notes.
"""
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd
import mpmath as mp

ZERO = 8

# ---------------------------------------------------------------- Q(sqrt2)
class K:
    __slots__ = ('p', 'q')
    def __init__(s, p=0, q=0): s.p = F(p); s.q = F(q)
    @staticmethod
    def c(x): return x if isinstance(x, K) else K(x)
    def __add__(s, o): o = K.c(o); return K(s.p + o.p, s.q + o.q)
    __radd__ = __add__
    def __neg__(s): return K(-s.p, -s.q)
    def __sub__(s, o): return s + (-K.c(o))
    def __rsub__(s, o): return K.c(o) - s
    def __mul__(s, o): o = K.c(o); return K(s.p*o.p + 2*s.q*o.q, s.p*o.q + s.q*o.p)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = K.c(o); n = o.p*o.p - 2*o.q*o.q
        assert n != 0
        r = s * K(o.p, -o.q); return K(r.p/n, r.q/n)
    def __rtruediv__(s, o): return K.c(o) / s
    def __eq__(s, o):
        if o is None: return False
        o = K.c(o); return s.p == o.p and s.q == o.q
    def __hash__(s): return hash((s.p, s.q))
    def sign(s):
        """exact sign of p + q sqrt2"""
        p, q = s.p, s.q
        if q == 0: return (p > 0) - (p < 0)
        if p == 0: return (q > 0) - (q < 0)
        if (p > 0) == (q > 0): return 1 if p > 0 else -1
        # opposite signs: compare p^2 and 2 q^2
        big_p = p*p > 2*q*q
        return (1 if p > 0 else -1) if big_p else (1 if q > 0 else -1)
    def __lt__(s, o): return (s - o).sign() < 0
    def __gt__(s, o): return (s - o).sign() > 0
    def __abs__(s): return -s if s.sign() < 0 else s
    def mpf(s): return mp.mpf(s.p.numerator)/s.p.denominator + mp.sqrt(2)*mp.mpf(s.q.numerator)/s.q.denominator
    def __repr__(s):
        return f'{s.p}' if s.q == 0 else (f'{s.q}*sqrt2' if s.p == 0 else f'({s.p}{"+" if s.q > 0 else ""}{s.q}*sqrt2)')

SQ2 = K(0, 1); A = K(1, 1); B = K(-1, 1)          # a = 1+sqrt2, b = sqrt2-1
INF = None
CONFIG = [K(0), K(1), K(-1), A, -A, B, -B, INF]

# ---------------------------------------------------------------- Cayley map
def cayley_exponent(x):
    """m with (x-i)/(x+i) = zeta^m, from the definition, checked numerically."""
    with mp.workdps(50):
        u = mp.mpc(1) if x is None else (x.mpf() - 1j) / (x.mpf() + 1j)
        m = int(mp.nint(mp.arg(u) / (mp.pi / 4))) % 8
        assert abs(u - mp.exp(2j * mp.pi * m / 8)) < mp.mpf(10)**-40, (x, u)
    return m
CAY = {x: cayley_exponent(x) for x in CONFIG}
assert sorted(CAY.values()) == list(range(8))

def pullback(p, q=None):
    """dlog((x-p)/(x-q)) as a form sum_m c_m dlog(u - zeta^m), u=(x-i)/(x+i).
    (x-p) = const*(u-u_p)/(u-1), so dlog(x-p) = dlog(u-u_p) - dlog(u-1); u_inf = 1."""
    f = Counter(); f[CAY[p]] += 1; f[CAY[q]] -= 1
    return {m: c for m, c in f.items() if c}

# ---------------------------------------------------------------- shuffle algebra
@lru_cache(maxsize=None)
def sh(u, v):
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = Counter()
    for w, c in sh(u[1:], v).items(): out[(u[0],) + w] += c
    for w, c in sh(u, v[1:]).items(): out[(v[0],) + w] += c
    return dict(out)

def vadd(*terms):
    out = {}
    for vec, c in terms:
        for w, a in vec.items():
            out[w] = out.get(w, 0) + c * a
    return {w: a for w, a in out.items() if a != 0}

def vmul(x, y):
    out = {}
    for u, a in x.items():
        for v, b in y.items():
            for w, c in sh(u, v).items():
                out[w] = out.get(w, 0) + a * b * c
    return {w: a for w, a in out.items() if a != 0}

def admissible(w): return bool(w) and w[0] != 0 and w[-1] != ZERO

# ---------------------------------------------------------------- path conversion
def piece_from_zero(forms, e):
    """int over the ray 0 -> zeta^e of forms f1..fk (f1 outermost) = sum of G(zeta^(m1-e),...;1)."""
    if not forms: return {(): 1}
    assert forms[0].get(e, 0) == 0, 'outermost form singular at the end point: would need regularisation'
    out = Counter()
    for choice in product(*[list(f.items()) for f in forms]):
        coeff = 1
        for m, c in choice: coeff *= c
        out[tuple((m - e) % 8 for m, c in choice)] += coeff
    return {w: c for w, c in out.items() if c}

def path_integral(forms, c, e):
    """Iterated integral of forms f1..fn (f1 outermost) along u: zeta^c -> 0 -> zeta^e.
    Chen: path = gamma1 (zeta^c -> 0) then gamma2 (0 -> zeta^e):
      I = sum_k I_gamma2(f1..fk) I_gamma1(f(k+1)..fn),  I_gamma1(g1..gm) = (-1)^m I_(0->zeta^c)(gm..g1).
    u = 0 is not a singular point of any form; no regularisation is used (asserted)."""
    n = len(forms); out = {}
    for k in range(n + 1):
        left = piece_from_zero(forms[:k], e)
        right = piece_from_zero(list(reversed(forms[k:])), c)
        out = vadd((out, 1), (vmul(left, right), (-1) ** (n - k)))
    assert all(admissible(w) for w in out), 'non-admissible word in conversion'
    return out

def strictly_between(x, lo, hi):
    lo, hi = (lo, hi) if lo < hi else (hi, lo)
    return x is not None and lo < x < hi

def convert_real_segment(xforms, x0, x1):
    """xforms: list of (p,q) meaning dlog((x-p)/(x-q)); integrate along the real segment x0 -> x1
    (finite).  Homotopy check: no singular point of any form lies strictly inside the segment
    (then the Cayley image arc and the radial path u(x0) -> 0 -> u(x1) are homotopic in the disc)."""
    assert x0 is not None and x1 is not None
    for p, q in xforms:
        for s in (p, q):
            assert not strictly_between(s, x0, x1), ('singular point inside the segment', s, x0, x1)
    return path_integral([pullback(p, q) for p, q in xforms], CAY[x0], CAY[x1])

# ---------------------------------------------------------------- cross ratios, atoms
def cross(x, p0, p1, pinf):
    """Moebius coordinate with p0 -> 0, p1 -> 1, pinf -> infinity (None = infinity)."""
    def d(u, v):
        if u is None and v is None: return K(0)
        if u is None: return K(1)
        if v is None: return K(-1)
        return u - v
    den = d(x, pinf) * d(p1, p0)
    return None if den == K(0) else d(x, p0) * d(p1, pinf) / den

def li_frames(z):
    """All (p0,p1,pinf,e) of configuration points with cross(e)=z, p0,e finite and the
    real segment [p0,e] free of pinf (and of p1; automatic for z<1)."""
    out = []
    for p0 in CONFIG:
        for p1 in CONFIG:
            for pinf in CONFIG:
                for e in CONFIG:
                    if len({repr(p0), repr(p1), repr(pinf), repr(e)}) < 4: continue
                    if p0 is None or e is None: continue
                    if cross(e, p0, p1, pinf) != z: continue
                    if strictly_between(pinf, p0, e) or strictly_between(p1, p0, e): continue
                    out.append((p0, p1, pinf, e))
    return out

def li_vector(n, z, frame):
    """Li_n(z) = -G(0^(n-1),1;z) in y = cross(x); dlog y = dlog((x-p0)/(x-pinf)),
    dlog(y-1) = dlog((x-p1)/(x-pinf)); path x: p0 -> e on the real line."""
    p0, p1, pinf, e = frame
    v = convert_real_segment([(p0, pinf)] * (n - 1) + [(p1, pinf)], p0, e)
    return {w: -c for w, c in v.items()}

LOG2 = {(4,): F(1)}                        # G(-1;1) = log 2
IPI = {(2,): F(2), (6,): F(-2)}            # 2(G(i;1)-G(-i;1)) = 2(log(1+i)-log(1-i)) = i pi
PI2 = {w: -c for w, c in vmul(IPI, IPI).items()}
ZETA3 = {(ZERO, ZERO, 0): F(-1)}           # zeta(3) = Li_3(1) = -G(0,0,1;1)
def LOGA(): return convert_real_segment([(K(0), None)], K(1), A)   # log a = int_1^a dx/x

# ---------------------------------------------------------------- H dictionary (comparison only)
def g_of_h(hw):
    sign = 1; out = []
    for j in hw:
        if j == -1: out.append(ZERO)
        else: out.append((-j) % 8); sign = -sign
    return tuple(out), sign
def h_of_g(gw):
    sign = 1; out = []
    for a in gw:
        if a == ZERO: out.append(-1)
        else: out.append((-a) % 8); sign = -sign
    return tuple(out), sign
def to_H(vec):
    out = {}
    for w, c in vec.items():
        hw, s = h_of_g(w); out[hw] = out.get(hw, 0) + s * c
    return {w: c for w, c in out.items() if c != 0}
def from_H(vec):
    out = {}
    for w, c in vec.items():
        gw, s = g_of_h(tuple(w)); out[gw] = out.get(gw, 0) + s * c
    return {w: c for w, c in out.items() if c != 0}

def primitive_H(vec):
    """Normalisation of the notes' rows (applied in H notation): integer coefficients with gcd 1,
    words sorted, first coefficient positive."""
    vec = {w: F(c) for w, c in vec.items() if c != 0}
    if not vec: return ()
    den = 1
    for c in vec.values(): den = den * c.denominator // gcd(den, c.denominator)
    items = sorted((w, int(c * den)) for w, c in vec.items())
    g = 0
    for w, c in items: g = gcd(g, c)
    if items[0][1] < 0: g = -g
    return tuple((w, c // g) for w, c in items)

# ---------------------------------------------------------------- numerics
def load_values(path, dps=140):
    import json
    mp.mp.dps = dps
    d = json.load(open(path))['values']
    return {tuple(int(x) for x in k.split()): mp.mpc(v[0], v[1]) for k, v in d.items()}
def mpc_of(c):
    if isinstance(c, K): return c.mpf()
    if isinstance(c, F): return mp.mpf(c.numerator) / c.denominator
    return mp.mpf(c)
def num(vec, vals):
    s = mp.mpc(0)
    for w, c in vec.items():
        s += mpc_of(c) * (1 if not w else vals[w])
    return s
