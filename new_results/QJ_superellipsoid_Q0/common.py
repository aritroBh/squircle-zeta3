# common.py -- shared exact tools for the QJ job.
#   * Q2: exact arithmetic in Q(sqrt2) (pairs of Fractions)
#   * logS(x): exact log of a positive S-unit of Q(sqrt2), S={2,3}, as (i/2)*lam + j*l3 + k*L
#   * sympy symbols lam=log2, l3=log3, L=log(1+sqrt2), al=arctan(1/sqrt2), Li2 = unevaluated function
#   * num(expr): numerical value (mpmath) of an expression in these symbols
from fractions import Fraction as Fr
import math
import sympy as sp
import mpmath as mp

lam, l3, L, al = sp.symbols('lambda ell3 L alpha_c', real=True)
PI = sp.pi
Li2 = sp.Function('Li2')          # unevaluated real dilogarithm Li2(x), x < 1 (x in Q(sqrt2))
R2 = sp.sqrt(2)


class Q2:
    __slots__ = ('c', 'd')       # c + d*sqrt2

    def __init__(self, c, d=0):
        self.c = Fr(c); self.d = Fr(d)

    def __add__(s, o): o = q2(o); return Q2(s.c + o.c, s.d + o.d)
    __radd__ = __add__
    def __neg__(s): return Q2(-s.c, -s.d)
    def __sub__(s, o): return s + (-q2(o))
    def __rsub__(s, o): return q2(o) - s
    def __mul__(s, o): o = q2(o); return Q2(s.c*o.c + 2*s.d*o.d, s.c*o.d + s.d*o.c)
    __rmul__ = __mul__
    def conj(s): return Q2(s.c, -s.d)
    def norm(s): return s.c*s.c - 2*s.d*s.d
    def inv(s):
        n = s.norm(); assert n != 0
        return Q2(s.c/n, -s.d/n)
    def __truediv__(s, o): return s*q2(o).inv()
    def __rtruediv__(s, o): return q2(o)*s.inv()
    def __pow__(s, k):
        r = Q2(1); bse = s if k >= 0 else s.inv()
        for _ in range(abs(k)): r = r*bse
        return r
    def __eq__(s, o): o = q2(o); return s.c == o.c and s.d == o.d
    def __hash__(s): return hash((s.c, s.d))
    def mp(s): return mp.mpf(s.c.numerator)/s.c.denominator + mp.mpf(s.d.numerator)/s.d.denominator*mp.sqrt(2)
    def __float__(s): return float(s.c) + float(s.d)*math.sqrt(2)
    def sign(s):
        v = float(s)
        if abs(v) > 1e-9: return 1 if v > 0 else -1
        return 1 if s.mp() > 0 else (-1 if s.mp() < 0 else 0)
    def __lt__(s, o): return (s - q2(o)).sign() < 0
    def __gt__(s, o): return (s - q2(o)).sign() > 0
    def __le__(s, o): return not s > o
    def __ge__(s, o): return not s < o
    def sym(s): return sp.Rational(s.c.numerator, s.c.denominator) + sp.Rational(s.d.numerator, s.d.denominator)*R2
    def __repr__(s):
        if s.d == 0: return str(s.c)
        if s.c == 0: return "%s*r2" % s.d
        return "(%s%s%s*r2)" % (s.c, '+' if s.d > 0 else '', s.d)


def q2(x):
    if isinstance(x, Q2): return x
    if isinstance(x, (int, Fr)): return Q2(x)
    # sympy number
    e = sp.expand(sp.radsimp(sp.sympify(x)))
    d = e.coeff(R2); c = sp.expand(e - d*R2)
    assert c.is_Rational and d.is_Rational, x
    return Q2(Fr(int(c.p), int(c.q)), Fr(int(d.p), int(d.q)))


A = Q2(1, 1); B = Q2(-1, 1); SQ2 = Q2(0, 1)


def vfac(n, p):
    n = abs(n); assert n != 0
    num, den = n.numerator, n.denominator; k = 0
    while num % p == 0: num //= p; k += 1
    while den % p == 0: den //= p; k -= 1
    return k, Fr(num, den)


def sunit(x):
    """x = sign * sqrt2^i * 3^j * a^k exactly; returns (sign,i,j,k) or None."""
    x = q2(x)
    if x == 0: return None
    N = x.norm()
    i, rest = vfac(N, 2); j2, rest = vfac(rest, 3)
    if rest != 1 or j2 % 2: return None
    j = j2//2
    y = x/(SQ2**i * Q2(3)**j)
    v = abs(float(y))
    k = round(math.log(v)/math.log(1 + math.sqrt(2)))
    for sg in (1, -1):
        if y == A**k*sg: return (sg, i, j, k)
    return None


def logS(x):
    """exact log of a positive S-unit x as a sympy expression in lam, l3, L"""
    t = sunit(x)
    assert t is not None and t[0] == 1, ("not a positive S-unit", x)
    _, i, j, k = t
    return sp.Rational(i, 2)*lam + j*l3 + k*L


def li2(x):
    return Li2(q2(x).sym())


NUMSUBS = None
def num(expr, dps=None):
    """numerical value of an expression in lam,l3,L,al,pi and Li2(x) with x < 1"""
    e = sp.sympify(expr)
    vals = {lam: mp.log(2), l3: mp.log(3), L: mp.log(1 + mp.sqrt(2)), al: mp.atan(1/mp.sqrt(2))}
    res = mp.mpf(0)
    e = sp.expand(e)
    for term in sp.Add.make_args(e):
        f = mp.mpf(1)
        for fac in sp.Mul.make_args(term):
            base, ex = fac.as_base_exp()
            if isinstance(base, sp.Function) and base.func == Li2:
                xv = q2(base.args[0]).mp(); assert xv < 1
                bv = mp.polylog(2, xv)
            elif base in vals: bv = vals[base]
            elif base == sp.pi: bv = mp.pi
            else: bv = mp.mpf(sp.N(base, mp.mp.dps + 10))
            f *= bv**int(ex) if ex.is_Integer else bv**mp.mpf(sp.N(ex, mp.mp.dps + 10))
        res += f
    return res


def show(name, val, ref, tol_digits=60):
    d = abs(val - ref)
    ok = d < mp.mpf(10)**(-tol_digits)
    print("  %-58s |diff| = %s %s" % (name, mp.nstr(d, 3), "OK" if ok else "FAIL"))
    assert ok, name
    return d
