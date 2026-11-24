"""Referee step 2b: own exact evaluation of T, W_3, W_3' (kernels re-derived in s2_ibp.py) and own
exact reduction of the resulting Li_2 / Li_3 values with classical functional equations.

Method (written independently of the project scripts):
  * split [b,a] at 1;  every term is int (l_x - k_x)(l_y - k_y) dv/(v-z), l_x = log|v-x|;
  * polarisation with the endpoint letter in the numerator of w=(v-x)/(v-y);
  * primitives of (kappa+log|s|)^n/(s-1), n=1,2, evaluated exactly at S-unit endpoints;
  * the pole of the kernel at the endpoint v=a is handled by shifting l_x -> l_x - log|a-x|.
Then the Li values are reduced by exact linear algebra over Q(sqrt2) using instances of
(R),(I),(La),(D),(A) at weight 2 and (I3),(T3),(D3),(K) at weight 3.
Usage: python s3_weight3.py [full|note]   ('note' = Kummer only at the four points of the note)
"""
import sys
import itertools
import sympy as sp
import mpmath as mp
from q2 import Q2, A, B, R2, ONE, ZERO, logS, lam, Lg

mp.mp.dps = 60
MODE = sys.argv[1] if len(sys.argv) > 1 else 'note'
pi2 = sp.pi**2
z3 = sp.Symbol('zeta3')

# ---------------- polylog symbols
REG = {}   # (n, Q2) -> symbol


def Li(n, u):
    if u == ZERO:
        return sp.Integer(0)
    key = (n, u)
    if key not in REG:
        REG[key] = sp.Symbol(f'Li{n}[{u}]')
    return REG[key]


def lg(u):          # log|u| exactly
    return logS(u)


# ---------------- primitives of log|s|/(s-1), log^2|s|/(s-1)
def region(smid):
    if smid < 0:
        return 'neg'
    return '01' if smid < 1 else 'gt1'


def P1(s, reg):
    if s == 'zero':
        return -pi2/6 if reg != 'neg' else 0
    if s == 'one':
        return 0
    if reg == 'neg':
        return lg(-s)*lg(1 - s) + Li(2, s)
    return -Li(2, 1 - s)


def P2(s, reg):
    if s == 'zero':
        return 0
    if s == 'one':
        return -2*z3
    if reg == '01':
        return lg(s)**2*lg(1 - s) + 2*lg(s)*Li(2, s) - 2*Li(3, s)
    if reg == 'neg':
        return lg(-s)**2*lg(1 - s) + 2*lg(-s)*Li(2, s) - 2*Li(3, s)
    r = ONE/s
    return -lg(r)**3/3 + lg(r)**2*lg(1 - r) + 2*lg(r)*Li(2, r) - 2*Li(3, r)


def Eprim(n, s, kappa, reg):
    """primitive of (kappa + log|s|)^n /(s-1), n = 1, 2, at s (exact) in region reg."""
    if s == 'one':
        assert sp.simplify(kappa) == 0, 'divergent endpoint'
        l1 = 0
    elif s == 'zero':
        l1 = 0
    else:
        l1 = lg(1 - s)
    if n == 1:
        return kappa*l1 + P1(s, reg)
    return kappa**2*l1 + 2*kappa*P1(s, reg) + P2(s, reg)


def tag(s):
    if s == ZERO:
        return 'zero'
    if s == ONE:
        return 'one'
    return s


def prim_between(n, sfun, kappa, p, q):
    mid = (p + q)/2
    reg = region(sfun(mid))
    return Eprim(n, tag(sfun(q)), kappa, reg) - Eprim(n, tag(sfun(p)), kappa, reg)


# ---------------- elementary integrals on I = [p, q]
def lx_at(x, v):
    return lg(v - x)


def int_const(z, p, q):
    assert z != p and z != q
    return lg(q - z) - lg(p - z)


def int_lin(x, kx, z, p, q):
    """int_p^q (l_x - kx) dv/(v-z)"""
    if x == z:
        assert x != p and x != q
        return ((lx_at(x, q) - kx)**2 - (lx_at(x, p) - kx)**2)/2
    sfun = lambda v: (v - x)/(z - x)
    return prim_between(1, sfun, lg(z - x) - kx, p, q)


def int_sq(x, kx, z, p, q):
    """int_p^q (l_x - kx)^2 dv/(v-z)"""
    if x == z:
        assert x != p and x != q
        return ((lx_at(x, q) - kx)**3 - (lx_at(x, p) - kx)**3)/3
    sfun = lambda v: (v - x)/(z - x)
    return prim_between(2, sfun, lg(z - x) - kx, p, q)


def int_quad(x, kx, y, ky, z, p, q):
    """int_p^q (l_x - kx)(l_y - ky) dv/(v-z)"""
    if x == y:
        return int_sq(x, kx, z, p, q)
    if y in (p, q):                 # endpoint letter goes in the numerator of w
        x, kx, y, ky = y, ky, x, kx
    alpha = ky - kx                 # (l_x-kx)-(l_y-ky) = log|w| + alpha
    w = lambda v: (v - x)/(v - y)
    if y == z:                      # dlog|v-z| = -dw/(w-1)
        Wt = -prim_between(2, w, alpha, p, q)
    elif x == z:                    # dlog|v-z| = dw/w - dw/(w-1)
        assert x not in (p, q)
        Wt = ((lg(w(q)) + alpha)**3 - (lg(w(p)) + alpha)**3)/3 - prim_between(2, w, alpha, p, q)
    else:
        cc = (x - z)/(y - z)
        Wt = prim_between(2, lambda v: w(v)/cc, alpha + lg(cc), p, q) - prim_between(2, w, alpha, p, q)
    return (int_sq(x, kx, z, p, q) + int_sq(y, ky, z, p, q) - Wt)/2


# ---------------- the integrands
LT = ({B: 1, -A: 1, ZERO: -1}, -2*lam)          # log t
LMT = ({A: 1, -B: 1, ZERO: -1}, -2*lam)         # log(1-t)
INTERVALS = [(B, ONE), (ONE, A)]


def shifted(form, z):
    """write form = sum c_x (l_x - log|z-x|) when form(z) = 0 (pole of the kernel at endpoint z)."""
    coeffs, C = form
    assert sp.expand(sum(cx*lg(z - x) for x, cx in coeffs.items()) + C) == 0
    return {x: (cx, lg(z - x)) for x, cx in coeffs.items()}


def plain(form):
    coeffs, C = form
    return {x: (cx, 0) for x, cx in coeffs.items()}, C


def integrate(term_list):
    """term_list: list of (forms, kernel) with forms a tuple of 1 or 2 log-linear forms and
    kernel dict {z: coefficient (sympy)} meaning sum k_z/(v-z)."""
    total = 0
    for forms, kern in term_list:
        for z, kz in kern.items():
            endpoint_pole = (z == A)
            for (p, q) in INTERVALS:
                if endpoint_pole:
                    fs = [(shifted(f, A), 0) for f in forms]
                else:
                    fs = [plain(f) for f in forms]
                if len(forms) == 1:
                    (cf, C), = fs
                    val = C*int_const(z, p, q) if C != 0 else 0
                    for x, (cx, kx) in cf.items():
                        val += cx*int_lin(x, kx, z, p, q)
                else:
                    (cf1, C1), (cf2, C2) = fs
                    val = C1*C2*int_const(z, p, q) if (C1 != 0 and C2 != 0) else 0
                    for x, (cx, kx) in cf1.items():
                        for y, (cy, ky) in cf2.items():
                            val += cx*cy*int_quad(x, kx, y, ky, z, p, q)
                    if C1 != 0:
                        for y, (cy, ky) in cf2.items():
                            val += C1*cy*int_lin(y, ky, z, p, q)
                    if C2 != 0:
                        for x, (cx, kx) in cf1.items():
                            val += C2*cx*int_lin(x, kx, z, p, q)
                total += kz*val
    return sp.expand(total)


s2 = sp.sqrt(2)
# kernels exactly as re-derived in s2_ibp.py (partial fractions)
T_terms = [((LT,), {ZERO: sp.Rational(1, 2), A: sp.Rational(1, 2), -A: -1, -B: -sp.Rational(1, 2)})]
W3_terms = [((LT, LT), {ZERO: (4 - 3*s2)/8, -A: -1}),
            ((LT, LMT), {ZERO: s2/8})]
W3p_terms = [((LT, LT), {ZERO: (1 - s2)/4, A: sp.Rational(1, 4), -A: -sp.Rational(1, 4), -B: -sp.Rational(1, 2)}),
             ((LT, LMT), {ZERO: (1 + s2)/4, -A: -sp.Rational(1, 2)})]


# ---------------- numerical evaluation of expressions
def numval(expr):
    sub = {lam: mp.log(2), Lg: mp.log(1 + mp.sqrt(2)), z3: mp.zeta(3)}
    for (n, u), sym in REG.items():
        sub[sym] = mp.polylog(n, u.mp()) if u.mp() <= 1 else mp.re(mp.polylog(n, u.mp()))
    f = sp.lambdify(list(sub.keys()), expr, 'mpmath', dummify=True)
    return f(*sub.values())


def vint(f):
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
    return mp.quad(f, [b, 1, a])


def l1(v):
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
    return mp.log((v - b)*(v + a)/(4*v))


def l2(v):
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1
    return mp.log((a - v)*(v + b)/(4*v))


if __name__ == '__main__':
    r2 = mp.sqrt(2); a_ = 1 + r2; b_ = r2 - 1
    Tx = integrate(T_terms)
    W3x = integrate(W3_terms)
    W3px = integrate(W3p_terms)
    T_q = vint(lambda v: l1(v)*(1/(2*v) + 1/(2*(v - a_)) - 1/(v + a_) - 1/(2*(v + b_))))
    W3_q = vint(lambda v: l1(v)**2*((4 - 3*r2)/(8*v) - 1/(v + a_)) + r2/8*l1(v)*l2(v)/v)
    W3p_q = vint(lambda v: l1(v)**2*((1 - r2)/(4*v) + 1/(4*(v - a_)) - 1/(4*(v + a_)) - 1/(2*(v + b_)))
                 + l1(v)*l2(v)*((1 + r2)/(4*v) - 1/(2*(v + a_))))
    for name, ex, qv in [('T', Tx, T_q), ('W3', W3x, W3_q), ("W3'", W3px, W3p_q)]:
        syms = sorted([str(s) for s in ex.free_symbols if str(s).startswith('Li')])
        print(f'{name}: raw exact expression uses {sum(1 for s in syms if s.startswith("Li3"))} Li3 and '
              f'{sum(1 for s in syms if s.startswith("Li2"))} Li2 symbols; |raw - quadrature| =',
              mp.nstr(abs(numval(ex) - qv), 3))
    print('Li2 args:', sorted({str(u) for (n, u) in REG if n == 2}))
    print('Li3 args:', sorted({str(u) for (n, u) in REG if n == 3}))
