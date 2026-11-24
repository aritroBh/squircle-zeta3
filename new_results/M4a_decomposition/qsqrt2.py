"""Exact arithmetic in Q(sqrt2) and exact logarithms of S-units.

Every number that occurs in the M1 computation is an element of Q(sqrt2), and every logarithm that
occurs is the log of an S-unit  +-2^(j/2) * a^m,  a = 1 + sqrt2.  So
    log|z| = (j/2) * log 2 + m * L,     L = log(1 + sqrt2),
exactly. No floating point is used anywhere in this module except to pick the integer m, which is
then checked exactly.
"""
from fractions import Fraction as Fr
import math


class QS2:
    """p + q*sqrt2 with rational p, q."""
    __slots__ = ("p", "q")

    def __init__(self, p=0, q=0):
        self.p = Fr(p)
        self.q = Fr(q)

    @staticmethod
    def c(x):
        return x if isinstance(x, QS2) else QS2(x, 0)

    def __add__(s, o):
        o = QS2.c(o); return QS2(s.p + o.p, s.q + o.q)
    __radd__ = __add__

    def __neg__(s):
        return QS2(-s.p, -s.q)

    def __sub__(s, o):
        return s + (-QS2.c(o))

    def __rsub__(s, o):
        return QS2.c(o) - s

    def __mul__(s, o):
        o = QS2.c(o); return QS2(s.p * o.p + 2 * s.q * o.q, s.p * o.q + s.q * o.p)
    __rmul__ = __mul__

    def conj(s):
        return QS2(s.p, -s.q)

    def norm(s):
        return s.p * s.p - 2 * s.q * s.q

    def __truediv__(s, o):
        o = QS2.c(o); n = o.norm()
        if n == 0:
            raise ZeroDivisionError
        r = s * o.conj(); return QS2(r.p / n, r.q / n)

    def __rtruediv__(s, o):
        return QS2.c(o) / s

    def __pow__(s, k):
        r = QS2(1)
        base = s if k >= 0 else QS2(1) / s
        for _ in range(abs(k)):
            r = r * base
        return r

    def __eq__(s, o):
        o = QS2.c(o); return s.p == o.p and s.q == o.q

    def __hash__(s):
        return hash((s.p, s.q))

    def sign(s):
        """Exact sign of p + q sqrt2."""
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
        if p > 0:   # p>0, q<0
            return 1 if p * p > 2 * q * q else -1
        return 1 if 2 * q * q > p * p else -1     # p<0, q>0

    def __lt__(s, o):
        return (s - o).sign() < 0

    def __le__(s, o):
        return (s - o).sign() <= 0

    def __gt__(s, o):
        return (s - o).sign() > 0

    def __ge__(s, o):
        return (s - o).sign() >= 0

    def __abs__(s):
        return s if s.sign() >= 0 else -s

    def __float__(s):
        return float(s.p) + float(s.q) * math.sqrt(2)

    def __repr__(s):
        if s.q == 0:
            return str(s.p)
        if s.p == 0:
            return f"{s.q}*sqrt2"
        return f"({s.p}{'+' if s.q > 0 else '-'}{abs(s.q)}*sqrt2)"

    def sympy(s):
        import sympy as sp
        return sp.Rational(s.p.numerator, s.p.denominator) + sp.Rational(s.q.numerator, s.q.denominator) * sp.sqrt(2)


SQ2 = QS2(0, 1)
A = QS2(1, 1)      # a = 1 + sqrt2  (fundamental unit)
B = QS2(-1, 1)     # b = sqrt2 - 1 = 1/a = tan(pi/8)
LOG_A = math.log(1 + math.sqrt(2))


def sunit_log(z):
    """For an S-unit z = +-2^(j/2) a^m return (j, m) with log|z| = (j/2) log2 + m L.  Exact check."""
    z = abs(QS2.c(z))
    n = abs(z.norm())                       # = 2^j
    j = n.numerator.bit_length() - 1 if n >= 1 else -((n.denominator.bit_length()) - 1)
    if Fr(2) ** j != n:
        raise ValueError(f"{z} is not an S-unit (norm {n})")
    # z / sqrt2^j is a unit u with |N(u)| = 1; u = a^m
    u = z / (SQ2 ** j)
    m = round(math.log(float(u)) / LOG_A)
    if u != A ** m:
        raise ValueError(f"{z}: unit part {u} is not a power of 1+sqrt2")
    return j, m


if __name__ == "__main__":
    for z in [B, A, SQ2, QS2(2) - SQ2, QS2(1, 0) / 2, QS2(3, -2), QS2(1) - SQ2 / 2, (QS2(2) - SQ2) / 4]:
        print(z, sunit_log(z))
