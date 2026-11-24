"""Exact arithmetic in Q(sqrt2) and exact logarithms of S-units +-2^(j/2) a^m, a = 1+sqrt2."""
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp


class Q2:
    __slots__ = ('p', 'q')

    def __init__(self, p, q=0):
        self.p = Fr(p)
        self.q = Fr(q)

    def __add__(s, o):
        o = c(o); return Q2(s.p + o.p, s.q + o.q)
    __radd__ = __add__

    def __sub__(s, o):
        o = c(o); return Q2(s.p - o.p, s.q - o.q)

    def __rsub__(s, o):
        return c(o) - s

    def __neg__(s):
        return Q2(-s.p, -s.q)

    def __mul__(s, o):
        o = c(o); return Q2(s.p*o.p + 2*s.q*o.q, s.p*o.q + s.q*o.p)
    __rmul__ = __mul__

    def conj(s):
        return Q2(s.p, -s.q)

    def norm(s):
        return s.p*s.p - 2*s.q*s.q

    def __truediv__(s, o):
        o = c(o); n = o.norm(); t = s*o.conj(); return Q2(t.p/n, t.q/n)

    def __rtruediv__(s, o):
        return c(o)/s

    def __eq__(s, o):
        if isinstance(o, str):
            return False
        o = c(o); return s.p == o.p and s.q == o.q

    def __hash__(s):
        return hash((s.p, s.q))

    def mp(s):
        return mp.mpf(s.p.numerator)/s.p.denominator + mp.mpf(s.q.numerator)/s.q.denominator*mp.sqrt(2)

    def sign(s):
        # exact sign of p + q sqrt2
        p, q = s.p, s.q
        if q == 0:
            return (p > 0) - (p < 0)
        if p == 0:
            return (q > 0) - (q < 0)
        if p > 0 and q > 0:
            return 1
        if p < 0 and q < 0:
            return -1
        # opposite signs: compare p^2 with 2 q^2
        d = p*p - 2*q*q
        return (1 if p > 0 else -1) if d > 0 else (1 if q > 0 else -1)

    def __lt__(s, o):
        return (s - c(o)).sign() < 0

    def __le__(s, o):
        return (s - c(o)).sign() <= 0

    def __gt__(s, o):
        return (s - c(o)).sign() > 0

    def __ge__(s, o):
        return (s - c(o)).sign() >= 0

    def __abs__(s):
        return -s if s.sign() < 0 else s

    def sym(s):
        return sp.Rational(s.p.numerator, s.p.denominator) + sp.Rational(s.q.numerator, s.q.denominator)*sp.sqrt(2)

    def __repr__(s):
        return f'({s.p}{"+" if s.q >= 0 else "-"}{abs(s.q)}r2)'


def c(o):
    return o if isinstance(o, Q2) else Q2(o)


R2 = Q2(0, 1)
A = Q2(1, 1)
B = Q2(-1, 1)
ONE = Q2(1)
ZERO = Q2(0)

lam, Lg = sp.symbols('lambda L')


def logS(u):
    """exact log|u| for an S-unit u = +-2^(j/2) a^m; returns sympy j/2*lambda + m*L."""
    u = abs(u)
    assert u != 0
    n = u.norm()              # = +-2^j
    n = abs(n)
    j = 0
    while n.numerator % 2 == 0 and n.numerator > 1:
        n /= 2; j += 1
    while n.denominator % 2 == 0:
        n *= 2; j -= 1
    assert n == 1, ('not an S-unit', u)
    w = u / pw(R2, j) if j >= 0 else u * pw(R2, -j)
    m = int(mp.nint(mp.log(w.mp())/mp.log(1 + mp.sqrt(2))))
    if m >= 0:
        assert w == pw(A, m), ('not S-unit', u)
    else:
        assert w == pw(B, -m), ('not S-unit', u)
    return sp.Rational(j, 2)*lam + m*Lg


def pw(x, n):
    r = ONE
    for _ in range(n):
        r = r*x
    return r


Q2.__pow__ = lambda s, n: pw(s, n)
