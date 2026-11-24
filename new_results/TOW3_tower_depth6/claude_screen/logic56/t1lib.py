"""t1lib.py -- referee's own exact Q(i,sqrt2) machinery for the pointwise T1 test (no tow2lib/tow3lib code in the
recursion itself; the closed forms are imported only as sympy expressions and compiled to exact evaluators)."""
import sys, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp

class GQ:
    """Gaussian rational a+bi"""
    __slots__ = ('a', 'b')
    def __init__(s, a=0, b=0):
        s.a = a if isinstance(a, Fr) else Fr(a); s.b = b if isinstance(b, Fr) else Fr(b)
    def __add__(s, o):
        if not isinstance(o, GQ): o = GQ(o)
        return GQ(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __neg__(s): return GQ(-s.a, -s.b)
    def __sub__(s, o):
        if not isinstance(o, GQ): o = GQ(o)
        return GQ(s.a-o.a, s.b-o.b)
    def __rsub__(s, o): return GQ(o) - s
    def __mul__(s, o):
        if not isinstance(o, GQ):
            o = Fr(o); return GQ(s.a*o, s.b*o)
        return GQ(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def __truediv__(s, o):
        if not isinstance(o, GQ):
            o = Fr(o); return GQ(s.a/o, s.b/o)
        n = o.a*o.a + o.b*o.b
        return s * GQ(o.a/n, -o.b/n)
    def conj(s): return GQ(s.a, -s.b)
    def iszero(s): return s.a == 0 and s.b == 0
    def __eq__(s, o):
        if not isinstance(o, GQ): o = GQ(o)
        return s.a == o.a and s.b == o.b
    def __repr__(s): return f'({s.a}+{s.b}i)'

I_ = GQ(0, 1)
def ipow(n):
    return [GQ(1), GQ(0, 1), GQ(-1), GQ(0, -1)][n % 4]

class QSI:
    """x + y*sqrt2 with x,y Gaussian rationals"""
    __slots__ = ('x', 'y')
    def __init__(s, x=None, y=None):
        s.x = x if x is not None else GQ(0); s.y = y if y is not None else GQ(0)
    def __add__(s, o): return QSI(s.x+o.x, s.y+o.y)
    def __sub__(s, o): return QSI(s.x-o.x, s.y-o.y)
    def __neg__(s): return QSI(-s.x, -s.y)
    def mul(s, g): return QSI(s.x*g, s.y*g)         # scalar in Q(i)
    def conj(s): return QSI(s.x.conj(), s.y.conj())
    def __eq__(s, o): return s.x == o.x and s.y == o.y
    def iszero(s): return s.x.iszero() and s.y.iszero()
    def __repr__(s): return f'[{s.x} + {s.y} r2]'
ZERO = QSI()

def fromdata(t):
    # TOWN1 layout: (Re A, Re B, Im A, Im B) for A + B sqrt2
    return QSI(GQ(Fr(t[0]), Fr(t[2])), GQ(Fr(t[1]), Fr(t[3])))

# ---------- compile sympy rational function in (d,m) over Q(i) into exact evaluator ----------
d_s, m_s = sp.symbols('d m')
def _poly_terms(expr):
    P = sp.Poly(sp.expand(expr), d_s, m_s)
    out = []
    for (i, j), ce in P.as_dict().items():
        ce = sp.sympify(ce)
        re_, im_ = ce.as_real_imag()
        re_ = sp.nsimplify(re_); im_ = sp.nsimplify(im_)
        assert re_.is_Rational and im_.is_Rational, ce
        out.append((i, j, Fr(int(re_.p), int(re_.q)), Fr(int(im_.p), int(im_.q))))
    return out
class RF:
    def __init__(self, expr):
        expr = sp.together(expr)
        num, den = sp.fraction(expr)
        self.num = _poly_terms(num); self.den = _poly_terms(den)
        self.cache = {}
    @staticmethod
    def _ev(terms, dd, mm, conj):
        re_ = Fr(0); im_ = Fr(0)
        for i, j, cr, ci in terms:
            w = (dd**i) * (mm**j)
            re_ += cr*w; im_ += (-ci if conj else ci)*w
        return GQ(re_, im_)
    def __call__(self, dd, mm, conj=False):
        key = (dd, mm, conj)
        r = self.cache.get(key)
        if r is None:
            r = self._ev(self.num, dd, mm, conj) / self._ev(self.den, dd, mm, conj)
            self.cache[key] = r
        return r

class Family:
    """e(delta,s)/e00 as exact evaluators; exprs: {(delta,s): sympy expr in d,m}"""
    def __init__(self, exprs):
        self.rf = {k: RF(v) for k, v in exprs.items()}
    def __call__(self, delta, s, dd, mm, conj=False):
        if s > delta//2: return GQ(0)
        return self.rf[(delta, s)](dd, mm, conj)

def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r
def e00_prime(p, q):
    """rho_pq(0) / (sqrt2/4)   [K_d (-1)^p (1+[p>q]) binom(d,p), K_d = (-1)^d (2d-3)!! i^{d+1} / (2 sqrt2 d!)]"""
    dd = p+q
    val = Fr((-1)**dd * dfact(2*dd-3), factorial(dd)) * (-1)**p * (2 if p > q else 1) * comb(dd, p)
    return ipow(dd+1) * val
def E00(p, q):
    return QSI(GQ(0), e00_prime(p, q) * Fr(1, 4))   # (sqrt2/4)*K'

_data = None
def load_data(path='<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl'):
    global _data
    if _data is None:
        _data = pickle.load(open(path, 'rb'))
    return _data
def data_get(p, q, k, delta):
    H = load_data()
    return fromdata(H[(p, q)][delta][k])
def data_has(p, q, k, delta):
    H = load_data()
    return (p, q) in H and delta in H[(p, q)] and k in H[(p, q)][delta]

# ---------- T1 (Proposition T1 of TOW1, re-implemented from its statement) ----------
class Taus:
    def __init__(self, t1, t2): self.t1 = t1; self.t2 = t2   # lists of GQ

def T1rhs(get, taus, p, q, k, delta, fold_sign=1, binom_fn=None):
    d = p+q
    tot = ZERO
    for dp in range(delta+1):
        c1 = taus.t1[delta-dp] * Fr(p+1, 2*d-dp-1)
        tot = tot + get(p+1, q, k, dp).mul(c1)
        if p > q:
            c2 = taus.t2[delta-dp] * Fr(q+1, 2*d-dp-1)
            tot = tot + get(p, q+1, k, dp).mul(c2)
    if q == p-1:
        mm = 4*p-4-delta
        for mu in range(mm, 4*p-4+1):
            g = ZERO
            for dp in range(0, 4*p-4-mu+1):
                g = g - get(p, p, k, dp).mul(taus.t1[4*p-4-mu-dp] / GQ(4*p-3-dp))
            bn = comb(mu-2, mu-mm) if binom_fn is None else binom_fn(mu-2, mu-mm)
            f = _mi_pow(mu) * ipow(mm-2) * bn          # (-i)^mu * i^(m-2) * binom(mu-2, mu-m)
            tot = tot + g.conj().mul(f * (-p) * fold_sign)
    return tot
def _mi_pow(n):
    # (-i)^n
    return [GQ(1), GQ(0, -1), GQ(-1), GQ(0, 1)][n % 4]
