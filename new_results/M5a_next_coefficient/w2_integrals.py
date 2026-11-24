"""Exact evaluation of  int_b^a Lambda(v) R(v) dv  where
   Lambda(v) = sum_x c_x log|v - x| + const   (x real S-units of Q(sqrt2)),
   R(v)      = rational function over Q(sqrt2) with only simple real poles,
as a W (weight-2 combination with formal Li2's), using the exact antiderivative of log|1-w|/w.
Poles at an endpoint (v = a or v = b) are allowed only where Lambda vanishes there; they are
handled by subtracting Lambda(endpoint) = 0 term by term (each subtracted integral converges).
"""
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
from qsqrt2 import QS2, SQ2, A, B
from dilog_exact import W, LOG, mul, PI2, LI2, float_exact

ONE = QS2(1)
LO, HI = B, A
v = sp.Symbol('v', positive=True)      # same symbol as in mk_integrands.py
S2 = sp.sqrt(2)


def to_qs2(expr):
    """exact sympy number in Q(sqrt2) -> QS2"""
    e = sp.nsimplify(sp.radsimp(sp.expand(expr)), [S2]) if False else sp.expand(sp.radsimp(expr))
    p = e.subs(S2, 0)
    q = sp.expand((e - p) / S2)
    assert q.free_symbols == set() and p.free_symbols == set()
    assert sp.simplify(e - (p + q * S2)) == 0
    p, q = sp.Rational(p), sp.Rational(q)
    return QS2(Fr(int(p.p), int(p.q)), Fr(int(q.p), int(q.q)))


def K(wv):
    """antiderivative of log|1-w|/w (continuous on each side of w = 0, K(0) = 0)."""
    if wv <= ONE:
        return -LI2(wv) if wv != QS2(0) else W()
    lw = LOG(wv)
    return mul(lw, lw).scale(Fr(1, 2)) + LI2(ONE / wv) - PI2.scale(Fr(1, 3))


def I(x, y):
    """int_b^a log|v-x| dv/(v-y),  y outside [b,a]."""
    assert not (LO <= y <= HI)
    if x == y:
        l1, l0 = LOG(HI - y), LOG(LO - y)
        return (mul(l1, l1) - mul(l0, l0)).scale(Fr(1, 2))
    w0, w1 = (LO - y) / (x - y), (HI - y) / (x - y)
    return mul(LOG(x - y), LOG(w1 / w0)) + K(w1) - K(w0)


def I_reg(x, y):
    """int_b^a (log|v-x| - log|y-x|) dv/(v-y),  y in {a, b} an endpoint, x != y."""
    assert y in (LO, HI) and x != y
    w0, w1 = (LO - y) / (x - y), (HI - y) / (x - y)
    return K(w1) - K(w0)


def partial_fractions(R):
    """R (sympy rational in v over Q(sqrt2), simple real poles) -> (poly, {pole: residue}) with QS2 entries."""
    terms = sp.Add.make_args(sp.expand(R))
    def _is_simple(term):
        num, den = sp.fraction(sp.together(term))
        return (not sp.expand(num).has(v)) and den.has(v) and sp.Poly(sp.expand(den), v).degree() == 1
    if all(_is_simple(tm) for tm in terms):
        pf = sp.expand(R)                      # already a sum of simple fractions
    else:
        pf = sp.apart(sp.together(R), v, extension=S2)
    poly = sp.Integer(0)
    res = {}
    for term in sp.Add.make_args(sp.expand(pf)):
        num, den = sp.fraction(sp.together(term))
        if not den.has(v):
            poly += term
            continue
        dpol = sp.Poly(sp.expand(den), v)
        assert dpol.degree() == 1, f"non-simple pole term {term}"
        lc = dpol.LC()
        y = to_qs2(-dpol.nth(0) / lc)
        k = to_qs2(num / lc)
        res[y] = res.get(y, QS2(0)) + k
    return sp.expand(poly), res


def lam_value_at(lam, y):
    """Lambda(y) as W (requires y not a log-letter of Lambda with nonzero coefficient)."""
    coeffs, const = lam
    out = const
    for x, c in coeffs.items():
        assert x != y
        out = out + LOG(y - x).scale(c)
    return out


def integrate(lam, R):
    """lam = ({x: c_x (QS2)}, const W);  R sympy rational.  Returns W = int_b^a Lambda R dv."""
    coeffs, const = lam
    poly, res = partial_fractions(R)
    assert poly == 0, f"polynomial part {poly} not supported"
    total = W()
    for y, r in res.items():
        if LO <= y <= HI:
            assert y in (LO, HI), y
            assert coeffs.get(y, QS2(0)) == QS2(0)
            val0 = lam_value_at(lam, y)
            assert val0.d == {}, f"Lambda does not vanish at endpoint pole {y}: {val0}"
            piece = W()
            for x, c in coeffs.items():
                piece = piece + I_reg(x, y).scale(c)
        else:
            piece = W()
            for x, c in coeffs.items():
                piece = piece + I(x, y).scale(c)
            piece = piece + mul(const, LOG((HI - y) / (LO - y)))
        total = total + piece.scale(r)
    return total


# the two logarithms of the problem, in the v-variable:
#   log t     = log(v-b) + log(v+a) - log v - 2 log 2
#   log(1-t)  = log(a-v) + log(v+b) - log v - 2 log 2
LOG_T = ({B: ONE, -A: ONE, QS2(0): -ONE}, LOG(QS2(4)).scale(-1))
LOG_1MT = ({A: ONE, -B: ONE, QS2(0): -ONE}, LOG(QS2(4)).scale(-1))
