"""qs3_symbol.py -- full number-level Goncharov symbol machinery over Q(sqrt3),
constant-retaining, exactly as A3_W4_proof/w4_symbol.py but with the Q(sqrt3) S-unit
factorisation. QS3 = (Fraction p, Fraction q) = p + q*sqrt(3).

V basis (log alphabet), slots:
  0: X = log(2+sqrt3)   (fundamental unit, norm +1)
  1: Y = log(1+sqrt3)   (norm -2)
  2: Z = log(sqrt3)     (norm -3)
  3: W = log(5+2*sqrt3) (norm 13; the prime 13 enters via the letter difference p1-p3)
S-unit group (S = primes above 2,3 and the prime (5+2sqrt3) above 13) mod torsion:
  <U=2+r3, V=1+r3, r3, P13=5+2r3>; |N| = 1,2,3,13. log 2 = 2Y-X, log 3 = 2Z.
"""
from fractions import Fraction as F
from functools import cache
from collections import defaultdict

class Q3V(tuple):
    """p + q*sqrt(3) as a value type: supports +, -, *, == (incl. with int/Fraction 0)."""
    def __new__(cls, p=0, q=0):
        if isinstance(p, (tuple, list)): p, q = p
        return super().__new__(cls, (F(p), F(q)))
    def __add__(self, o):
        if isinstance(o, (int, F)): o = Q3V(o)
        return Q3V(self[0]+o[0], self[1]+o[1])
    __radd__ = __add__
    def __sub__(self, o):
        if isinstance(o, (int, F)): o = Q3V(o)
        return Q3V(self[0]-o[0], self[1]-o[1])
    def __rsub__(self, o): return Q3V(o) - self
    def __neg__(self): return Q3V(-self[0], -self[1])
    def __mul__(self, o):
        if isinstance(o, (int, F)): return Q3V(self[0]*o, self[1]*o)
        return Q3V(self[0]*o[0]+3*self[1]*o[1], self[0]*o[1]+self[1]*o[0])
    __rmul__ = __mul__
    def __eq__(self, o):
        if isinstance(o, (int, F)) and o == 0:
            return self[0] == 0 and self[1] == 0
        return tuple.__eq__(self, o)
    def __hash__(self): return tuple.__hash__(self)

ZERO = Q3V(0); ONE = Q3V(1); MONE = Q3V(-1)
def q3(r, s=0): return Q3V(F(r), F(s))
def q3add(a, b): return Q3V(a[0]+b[0], a[1]+b[1])
def q3sub(a, b): return Q3V(a[0]-b[0], a[1]-b[1])
def q3neg(a): return Q3V(-a[0], -a[1])
def q3mul(a, b): return Q3V(a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def q3div(a, b):
    n = b[0]*b[0]-3*b[1]*b[1]
    return Q3V((a[0]*b[0]-3*a[1]*b[1])/n, (a[1]*b[0]-a[0]*b[1])/n)
def q3pow(a, n):
    r = ONE
    if n < 0: a = q3div(ONE, a); n = -n
    for _ in range(n): r = q3mul(r, a)
    return r
def q3norm(a): return a[0]*a[0]-3*a[1]*a[1]
def q3sign(a):
    """sign of p+q*sqrt(3) as a real number; exact."""
    p, q = a
    if q == 0: return (p > 0) - (p < 0)
    if p == 0: return (q > 0) - (q < 0)
    if (p > 0) == (q > 0): return 1 if p > 0 else -1
    c = p*p - 3*q*q
    if c == 0: raise ValueError("sqrt3 rational")
    return (1 if p > 0 else -1) if c > 0 else ((1 if q > 0 else -1))
def q3abs(a): return a if q3sign(a) > 0 else q3neg(a)
def q3conj(a): return (a[0], -a[1])

U  = (F(2), F(1))        # 2+sqrt3
V  = (F(1), F(1))        # 1+sqrt3
R3 = (F(0), F(1))        # sqrt3
P13 = (F(5), F(2))       # 5+2sqrt3, norm 13

@cache
def logvec(z):
    """log|z| = a*X + b*Y + c*Z + d*W, exact; raises if z not in the S-unit group."""
    z = q3abs(z)
    if z == ZERO:
        return (F(0), F(0), F(0), F(0))   # unit-tangent constant term (coincident points)
    nrm = abs(q3norm(z))
    p, q = nrm.numerator, nrm.denominator
    def strip(n, p_):
        e = 0
        while n % p_ == 0: n //= p_; e += 1
        return n, e
    p, e2 = strip(p, 2); p, e3 = strip(p, 3); p, e13 = strip(p, 13)
    q, f2 = strip(q, 2); q, f3 = strip(q, 3); q, f13 = strip(q, 13)
    if p != 1 or q != 1:
        raise ValueError(f"not an S-unit: {z}, norm={nrm}")
    b = e2 - f2; c = e3 - f3; d = e13 - f13
    u = q3div(z, q3mul(q3mul(q3pow(V, b), q3pow(R3, c)), q3pow(P13, d)))
    # u must now be +/- U^a
    a = 0
    neg = q3sign(u) < 0
    if neg: u = q3neg(u)
    while q3sign(q3sub(u, ONE)) > 0:
        u = q3div(u, U); a += 1
    while q3sign(q3sub(u, ONE)) < 0:
        u = q3mul(u, U); a -= 1
        if q3sign(q3sub(u, ONE)) > 0:
            raise ValueError(f"not an S-unit (unit part): {z}")
    assert u == ONE, f"not an S-unit (residual): {z}"
    # verify
    assert q3abs(q3mul(q3mul(q3pow(U, a), q3pow(V, b)), q3mul(q3pow(R3, c), q3pow(P13, d)))) == z
    return (F(a), F(b), F(c), F(d))

# ---------- symbol dict algebra ----------
def clean(d): return {k: v for k, v in d.items() if v != ZERO}
def add(*ts):
    out = defaultdict(lambda: ZERO)
    for t in ts:
        for k, v in t.items(): out[k] = q3add(out[k], v)
    return clean(out)
def scale(t, c): return clean({k: q3mul(k2, c) for k, k2 in ((k, v) for k, v in t.items())})
@cache
def shuffles(u, v):
    if not u: return (v,)
    if not v: return (u,)
    return tuple((u[0],)+w for w in shuffles(u[1:], v)) + tuple((v[0],)+w for w in shuffles(u, v[1:]))
def shuffle(s, t):
    out = defaultdict(lambda: ZERO)
    for u, c in s.items():
        for v, d in t.items():
            for w in shuffles(u, v): out[w] = q3add(out[w], q3mul(c, d))
    return clean(out)
def concat(s, t):
    out = defaultdict(lambda: ZERO)
    for u, c in s.items():
        for v, d in t.items(): out[u+v] = q3add(out[u+v], q3mul(c, d))
    return clean(out)
def log_symbol(z):
    return {(i,): q3(v) for i, v in enumerate(logvec(z)) if v}

@cache
def g_symbol(c, w, e):
    """full number-level Goncharov symbol of G_c(w; e), constant-retaining.
    c, e, letters of w: QS3. Word letters of the symbol: V slots."""
    if not w: return {(): ONE}
    if c == e: return {}
    out = {}
    points = (e,) + w + (c,)
    for i, x in enumerate(w, 1):
        ell = add(log_symbol(q3sub(x, points[i-1])), scale(log_symbol(q3sub(x, points[i+1])), MONE))
        out = add(out, concat(g_symbol(c, w[:i-1]+w[i:], e), ell))
    return out

def li_symbol(n, z):
    out = scale(log_symbol(q3sub(ONE, z)), MONE)
    for _ in range(n-1):
        out = concat(out, log_symbol(z))
    return out

def const_pow(slotvec, k):
    """symbol of (ell)^k where ell = sum slotvec[i]*V_i: k! * constant word on the vector..."""
    # returns shuffle-power of the weight-1 symbol
    import math
    base = {(i,): q3(c) for i, c in enumerate(slotvec) if c}
    out = {(): ONE}
    for _ in range(k): out = shuffle(out, base)
    return out
