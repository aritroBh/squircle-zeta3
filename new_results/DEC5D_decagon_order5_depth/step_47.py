
lib = '''"""dec5d_lib.py — exact quadratic-field arithmetic and generic Goncharov-symbol machinery
for the DEC5D depth test (regular decagon, Q(sqrt5)) with the octagon (Q(sqrt2)) as control.

Design mirrors the validated chain qsqrt2.py -> w4_symbol.py -> oct5d2_engine/symengine.py,
but is field-generic: a "kit" provides the exact number class, an exact S-unit logvec, the
G-letters and the log-alphabet slot names. The full number-level symbol keeps EVERY constant
entry (a_i - a_j), exactly as w4_symbol.g_symbol (the OCT5D2 lesson).

Obstruction: the antisymmetrised (3,2) cobracket component of a weight-5 symbol s,
  obs_{p,q}(s) = sum_{u in p, v in q} c_u d_v (s(u+v) - s(v+u)),
for p in a Hall basis of Lie_3(V), q in a Hall basis of Lie_2(V)  (d = dim V >= 2).
This is <delta_{3,2}([s]) - swap(delta_{2,3}([s])), p x q>: the projection modulo products.
For d = 2 it is exactly M4c's two pairings. The Hall pairing matrix (Lie elements as
functionals against Hall words) is checked invertible exactly, so all-zero <=> vanishing of
the (3,2) cobracket component in Lie_3(V) x Lie_2(V).
"""
from fractions import Fraction as F
from functools import cache
import math

# ---------------- exact Q(sqrt5) ----------------
class QS5:
    """p + q*sqrt5 with rational p, q."""
    __slots__ = ("p", "q")
    def __init__(self, p=0, q=0):
        self.p = F(p); self.q = F(q)
    @staticmethod
    def c(x):
        return x if isinstance(x, QS5) else QS5(x, 0)
    def __add__(s, o):
        o = QS5.c(o); return QS5(s.p + o.p, s.q + o.q)
    __radd__ = __add__
    def __neg__(s): return QS5(-s.p, -s.q)
    def __sub__(s, o): return s + (-QS5.c(o))
    def __rsub__(s, o): return QS5.c(o) - s
    def __mul__(s, o):
        o = QS5.c(o); return QS5(s.p*o.p + 5*s.q*o.q, s.p*o.q + s.q*o.p)
    __rmul__ = __mul__
    def conj(s): return QS5(s.p, -s.q)
    def norm(s): return s.p*s.p - 5*s.q*s.q
    def __truediv__(s, o):
        o = QS5.c(o); n = o.norm()
        if n == 0: raise ZeroDivisionError
        r = s*o.conj(); return QS5(r.p/n, r.q/n)
    def __rtruediv__(s, o): return QS5.c(o) / s
    def __pow__(s, k):
        r = QS5(1)
        base = s if k >= 0 else QS5(1)/s
        for _ in range(abs(k)): r = r*base
        return r
    def __eq__(s, o):
        o = QS5.c(o); return s.p == o.p and s.q == o.q
    def __hash__(s): return hash((s.p, s.q))
    def sign(s):
        p, q = s.p, s.q
        if q == 0: return (p > 0) - (p < 0)
        if p == 0: return (q > 0) - (q < 0)
        if p > 0 and q > 0: return 1
        if p < 0 and q < 0: return -1
        if p > 0: return 1 if p*p > 5*q*q else -1
        return 1 if 5*q*q > p*p else -1
    def __lt__(s, o): return (s - o).sign() < 0
    def __le__(s, o): return (s - o).sign() <= 0
    def __gt__(s, o): return (s - o).sign() > 0
    def __ge__(s, o): return (s - o).sign() >= 0
    def __abs__(s): return s if s.sign() >= 0 else -s
    def __float__(s): return float(s.p) + float(s.q)*math.sqrt(5)
    def __repr__(s):
        if s.q == 0: return str(s.p)
        if s.p == 0: return f"{s.q}*sqrt5"
        return f"({s.p}{'+' if s.q > 0 else '-'}{abs(s.q)}*sqrt5)"
    def sympy(s):
        import sympy as sp
        return sp.Rational(s.p.numerator, s.p.denominator) + sp.Rational(s.q.numerator, s.q.denominator)*sp.sqrt(5)

SQ5 = QS5(0, 1)
PHI = QS5(F(1,2), F(1,2))      # (1+sqrt5)/2
LOG_PHI = math.log((1+math.sqrt(5))/2)

def _v_p(n, p):
    n = abs(n); e = 0
    while n % p == 0: n //= p; e += 1
    return e, n

@cache
def logvec5(z):
    """log|z| = j*log2 + k*log(sqrt5) + m*log(phi), found and verified without floats.
    Returns (j, k, m) as Fractions. Slots: (log2, logsqrt5, logphi)."""
    z = abs(QS5.c(z))
    if z == QS5(0):
        return (F(0), F(0), F(0))     # unit-tangent constant term, not an ordinary log
    norm = abs(z.norm())
    a2, r = _v_p(norm.numerator, 2); a2d, rd = _v_p(norm.denominator, 2)
    j2n = a2 - a2d
    if j2n % 2: raise ValueError(f"odd 2-adic norm valuation: {z}")
    b5, r = _v_p(r, 5); b5d, rd = _v_p(rd, 5)
    k5 = b5 - b5d
    if r != 1 or rd != 1:
        raise ValueError(f"not an S-unit: {z}, norm={norm}")
    u = z / (QS5(2)**(j2n//2)) / (SQ5**k5)
    m = 0
    while u > QS5(1):
        u = u / PHI; m += 1
        if m > 400: raise ValueError(f"not an S-unit: {z}")
    while u < QS5(1):
        u = u * PHI; m -= 1
        if m < -400: raise ValueError(f"not an S-unit: {z}")
    assert u == QS5(1), f"unit part of {z} is not a power of phi"
    return (F(j2n//2), F(k5), F(m))

# ---------------- generic exact quadratic kit for Q(sqrt2) (octagon control) ----------------
import sys as _sys
_sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2, SQ2, A as A2
LOG_A2 = math.log(1+math.sqrt(2))

@cache
def logvec2(z):
    """log|z| = (j/2)*log2 + m*log(1+sqrt2): identical convention to w4_symbol.logvec,
    reimplemented here (cross-checked against it)."""
    z = abs(QS2.c(z))
    if z == QS2(0):
        return (F(0), F(0))
    norm = abs(z.norm())
    a2, r = _v_p(norm.numerator, 2); a2d, rd = _v_p(norm.denominator, 2)
    j = a2 - a2d
    if r != 1 or rd != 1:
        raise ValueError(f"not an S-unit: {z}, norm={norm}")
    u = z / (SQ2**j)
    m = 0
    while u > QS2(1):
        u = u / A2; m += 1
        if m > 400: raise ValueError(f"not an S-unit: {z}")
    while u < QS2(1):
        u = u * A2; m -= 1
        if m < -400: raise ValueError(f"not an S-unit: {z}")
    assert u == QS2(1), f"unit part of {z} not a power of 1+sqrt2"
    return (F(j, 2), F(m))

# ---------------- generic symbol machinery (mirrors w4_symbol.py) ----------------
def make_kit(Q, logvec, slots):
    ZERO, ONE = Q(0), Q(1)
    def clean(d): return {k: v for k, v in d.items() if v != ZERO}
    def add(*ts):
        out = {}
        for t in ts:
            for k, v in t.items(): out[k] = out.get(k, ZERO) + v
        return clean(out)
    def scale(t, c): return clean({k: v*c for k, v in t.items()})
    @cache
    def shuffles(u, v):
        if not u: return (v,)
        if not v: return (u,)
        return tuple((u[0],)+w for w in shuffles(u[1:], v)) + tuple((v[0],)+w for w in shuffles(u, v[1:]))
    def shuffle(s, t):
        out = {}
        for u, c in s.items():
            for v, d in t.items():
                for w in shuffles(u, v): out[w] = out.get(w, ZERO) + c*d
        return clean(out)
    def concat(s, t):
        out = {}
        for u, c in s.items():
            for v, d in t.items(): out[u+v] = out.get(u+v, ZERO) + c*d
        return clean(out)
    def log_symbol(z):
        return {(i,): Q(v) for i, v in enumerate(logvec(z)) if v}
    @cache
    def g_symbol(c, w, e):
        if not w: return {(): ONE}
        if c == e: return {}
        out = {}
        points = (e,) + w + (c,)
        for i, x in enumerate(w, 1):
            ell = add(log_symbol(x - points[i-1]), scale(log_symbol(x - points[i+1]), -1))
            out = add(out, concat(g_symbol(c, w[:i-1]+w[i:], e), ell))
        return out
    def li_symbol(n, z):
        out = scale(log_symbol(ONE - z), -1)
        for _ in range(n-1): out = concat(out, log_symbol(z))
        return out
    def const_pow(slot, k):
        return {(slot,)*k: Q(math.factorial(k))}
    def bracket(x, y): return add(concat(x, y), scale(concat(y, x), -1))
    return dict(Q=Q, ZERO=ZERO, ONE=ONE, logvec=logvec, slots=slots, clean=clean, add=add,
                scale=scale, shuffles=shuffles, shuffle=shuffle, concat=concat,
                log_symbol=log_symbol, g_symbol=g_symbol, li_symbol=li_symbol,
                const_pow=const_pow, bracket=bracket)

def hall_bases(d):
    """Hall basis of Lie_2 and Lie_3 over d letters (slot indices 0..d-1), as raw
    (list of (word, coeff)) pairs; words are tuples. Lie_3: [x_i,[x_j,x_k]], j<k, i>=j."""
    L2 = []
    for j in range(d):
        for k in range(j+1, d):
            L2.append([((j,k), 1), ((k,j), -1)])
    L3 = []
    for j in range(d):
        for k in range(j+1, d):
            for i in range(j, d):
                # [x_i, [x_j,x_k]] = x_i x_j x_k - x_i x_k x_j - x_j x_k x_i + x_k x_j x_i
                L3.append([((i,j,k),1), ((i,k,j),-1), ((j,k,i),-1), ((k,j,i),1)])
    return L2, L3

def make_obstruction(kit):
    Q, ZERO = kit['Q'], kit['ZERO']
    d = len(kit['slots'])
    L2r, L3r = hall_bases(d)
    conv = lambda terms: {w: Q(c) for w, c in terms}
    L2 = [conv(t) for t in L2r]; L3 = [conv(t) for t in L3r]
    def obstruction(s):
        return [[sum((c*dd*(s.get(u+v, ZERO) - s.get(v+u, ZERO))
                      for u, c in p.items() for v, dd in q.items()), ZERO)
                 for q in L2] for p in L3]
    return obstruction, L2, L3

def hall_duality_check(kit):
    """exact check: the pairing matrix <B(p), w'> (Lie_3 Hall elements as functionals vs
    Hall words) is invertible, so the obstruction pairings separate Lie_3 x Lie_2."""
    from fractions import Fraction as Fr
    d = len(kit['slots'])
    L2r, L3r = hall_bases(d)
    def wordset(terms): return {w for w, c in terms}
    # matrix for Lie_3 elements against their own words
    import itertools
    def det(mat):
        n = len(mat); M = [row[:] for row in mat]; dd = Fr(1)
        for i in range(n):
            piv = next((r for r in range(i, n) if M[r][i]), None)
            if piv is None: return Fr(0)
            if piv != i: M[i], M[piv] = M[piv], M[i]; dd = -dd
            dd *= M[i][i]
            for r in range(i+1, n):
                f = M[r][i]/M[i][i]
                for c in range(i, n): M[r][c] -= f*M[i][c]
        return dd
    M3 = [[Fr(sum(c for w, c in p if w == w2)) for w2, _ in [(w, 1) for p2 in L3r for w, c2 in p2[:1]]] for p in L3r]
    # simpler: pair each Lie element against the first word of each Hall element (its Hall word)
    hw3 = [p[0][0] for p in L3r]   # (i,j,k) is NOT the Hall word; Hall word of [x_i,[x_j,x_k]] is j k i
    hw3 = [(p[2][0]) for p in L3r] # the word (j,k,i): the actual Hall/Lyndon word
    M3 = [[Fr(dict(p).get(hw, 0)) for hw in hw3] for p in L3r]
    hw2 = [(p[0][0]) for p in L2r] # (j,k) is the Hall word for [x_j,x_k]
    M2 = [[Fr(dict(p).get(hw, 0)) for hw in hw2] for p in L2r]
    return det(M3), det(M2)

# ---------------- generic exact rewrite engine (mirrors oct5d2_engine.py) ----------------
def make_rewrite(Q, P_LET, Q_LET):
    ZERO, ONE = Q(0), Q(1)
    @cache
    def shuf(u, v):
        if not u: return {v: ONE}
        if not v: return {u: ONE}
        out = {}
        for w, c in shuf(u[1:], v).items(): out[(u[0],)+w] = out.get((u[0],)+w, ZERO) + c
        for w, c in shuf(u, v[1:]).items(): out[(v[0],)+w] = out.get((v[0],)+w, ZERO) + c
        return out
    def shuffle_words(parts):
        out = {(): ONE}
        for part in parts:
            nxt = {}
            for u, cu in out.items():
                for w, c in shuf(u, part).items():
                    nxt[w] = nxt.get(w, ZERO) + cu*c
            out = nxt
        return out
    @cache
    def IG(n, w):
        if not w: return {(): ONE/Q(n+1)}
        w1, wr = w[0], w[1:]
        out = {}
        c = (ONE - w1**(n+1)) / Q(n+1)
        if c != ZERO: out[w] = c
        coeffs = [w1**(n-j) for j in range(n+1)]
        for j in range(n+1):
            cj = coeffs[j]
            if cj == ZERO: continue
            for v, cv in IG(j, wr).items():
                out[v] = out.get(v, ZERO) - cj*cv/Q(n+1)
        return {k: v for k, v in out.items() if v != ZERO}
    def monomial_integral(n, a, b, c, d):
        parts = []
        if a: parts.append((ZERO,)*a)
        if b: parts.append((ONE,)*b)
        if c: parts.append((P_LET,)*c)
        if d: parts.append((Q_LET,)*d)
        sw = shuffle_words(parts)
        fac = Q(math.factorial(a)*math.factorial(b)*math.factorial(c)*math.factorial(d))
        out = {}
        for w, mult in sw.items():
            for v, cv in IG(n, w).items():
                out[v] = out.get(v, ZERO) + fac*mult*cv
        return {k: v for k, v in out.items() if v != ZERO}
    return dict(shuf=shuf, shuffle_words=shuffle_words, IG=IG, monomial_integral=monomial_integral)
'''
open("dec5d_lib.py","w").write(lib)
print("written dec5d_lib.py", len(lib))
import dec5d_lib as D
print("import OK; QS5 + PHI:", D.PHI, "=", float(D.PHI))
print("logvec5(phi) =", D.logvec5(D.PHI))
print("logvec5(2) =", D.logvec5(D.QS5(2)))
print("logvec5(sqrt5) =", D.logvec5(D.SQ5))
