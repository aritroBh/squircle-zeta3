"""cwfield.py -- scalar fields for the corner-family reducer.
make_quad(s): Q(theta), theta^2 = s (Fraction, non-square), elements a + b theta, a,b in Q.
QSI: Q(sqrt2, i) 4-tuple (copied interface from WTS fastred) for the squircle control c=1/sqrt2.
Common interface: + - * / ** ==, is_zero, is_real, conj (theta -> -theta; for QSI: i -> -i), scale(Fr).
"""
from fractions import Fraction as Fr

def make_quad(s):
    s = Fr(s)
    class Quad:
        __slots__ = ("a", "b")
        modulus = s
        def __init__(self, a=0, b=0):
            self.a = a if isinstance(a, Fr) else Fr(a)
            self.b = b if isinstance(b, Fr) else Fr(b)
        @staticmethod
        def c(x):
            return x if isinstance(x, Quad) else Quad(x)
        def __add__(self, o):
            o = Quad.c(o); return Quad(self.a + o.a, self.b + o.b)
        __radd__ = __add__
        def __neg__(self): return Quad(-self.a, -self.b)
        def __sub__(self, o): return self + (-Quad.c(o))
        def __rsub__(self, o): return Quad.c(o) + (-self)
        def scale(self, q): return Quad(self.a*Fr(q), self.b*Fr(q))
        def __mul__(self, o):
            o = Quad.c(o)
            return Quad(self.a*o.a + s*self.b*o.b, self.a*o.b + self.b*o.a)
        __rmul__ = __mul__
        def conj(self): return Quad(self.a, -self.b)
        def inv(self):
            n = self.a*self.a - s*self.b*self.b
            assert n != 0, "inv of 0"
            return Quad(self.a/n, -self.b/n)
        def __truediv__(self, o): return self * Quad.c(o).inv()
        def __rtruediv__(self, o): return Quad.c(o) * self.inv()
        def __pow__(self, k):
            r = Quad(1); base = self if k >= 0 else self.inv()
            for _ in range(abs(k)): r = r*base
            return r
        def is_zero(self): return self.a == 0 and self.b == 0
        def is_real(self): return self.b == 0
        def __eq__(self, o):
            o = Quad.c(o); return self.a == o.a and self.b == o.b
        def __hash__(self): return hash((self.a, self.b))
        def __repr__(self): return f"Q({self.a},{self.b})"
    return Quad


class QSI:
    """Q(sqrt2, i): p0 + p1 sqrt2 + p2 i + p3 sqrt2 i.  conj = i -> -i (theta-conjugation)."""
    __slots__ = ("p0", "p1", "p2", "p3")
    def __init__(s, p0=0, p1=0, p2=0, p3=0):
        s.p0 = Fr(p0); s.p1 = Fr(p1); s.p2 = Fr(p2); s.p3 = Fr(p3)
    @staticmethod
    def c(x):
        return x if isinstance(x, QSI) else QSI(x)
    def __add__(s, o):
        o = QSI.c(o); return QSI(s.p0+o.p0, s.p1+o.p1, s.p2+o.p2, s.p3+o.p3)
    __radd__ = __add__
    def __neg__(s):
        return QSI(-s.p0, -s.p1, -s.p2, -s.p3)
    def __sub__(s, o):
        return s + (-QSI.c(o))
    def __rsub__(s, o):
        return QSI.c(o) - s
    def scale(s, q):
        q = Fr(q); return QSI(s.p0*q, s.p1*q, s.p2*q, s.p3*q)
    def __mul__(s, o):
        o = QSI.c(o)
        a,b,c,d = s.p0,s.p1,s.p2,s.p3
        e,f,g,h = o.p0,o.p1,o.p2,o.p3
        return QSI(a*e + 2*b*f - c*g - 2*d*h,
                   a*f + b*e - c*h - d*g,
                   a*g + c*e + 2*b*h + 2*d*f,
                   a*h + d*e + b*g + c*f)
    __rmul__ = __mul__
    def conj(s):   # i -> -i
        return QSI(s.p0, s.p1, -s.p2, -s.p3)
    def inv(s):
        cs = s.conj()
        N = s * cs                     # in Q(sqrt2): p2=p3=0
        assert N.p2 == 0 and N.p3 == 0
        # invert N = p0 + p1 sqrt2
        n2 = N.p0*N.p0 - 2*N.p1*N.p1
        assert n2 != 0
        Ni = QSI(N.p0/n2, -N.p1/n2, 0, 0)
        return cs * Ni
    def __truediv__(s, o):
        return s * QSI.c(o).inv()
    def __rtruediv__(s, o):
        return QSI.c(o) / s
    def __pow__(s, k):
        r = QSI(1); base = s if k >= 0 else s.inv()
        for _ in range(abs(k)):
            r = r * base
        return r
    def is_zero(s):
        return s.p0 == 0 and s.p1 == 0 and s.p2 == 0 and s.p3 == 0
    def is_real(s):
        return s.p2 == 0 and s.p3 == 0
    def __eq__(s, o):
        o = QSI.c(o)
        return s.p0==o.p0 and s.p1==o.p1 and s.p2==o.p2 and s.p3==o.p3
    def __hash__(s):
        return hash((s.p0,s.p1,s.p2,s.p3))
    def __repr__(s):
        return f"QSI({s.p0},{s.p1},{s.p2},{s.p3})"
