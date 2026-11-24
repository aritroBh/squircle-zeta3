common = r'''
"""cor3_common.py -- shared infrastructure for the COR3k proof.

Design rule (per the brief): sympy is NEVER asked to simplify a polylogarithm
identity.  All reductions are done by explicit, hand-coded rules:
  * d/dx polylog(s, a(x)) = polylog(s-1, a) * a'/a   (sympy does this)
  * polylog(1, a)  ->  -log(1-a), with -log(1-a) expanded BY HAND into the
    basis  {log x, log(1-x), log(1+x), log 2}  (all factors positive on (0,1)).
  * vanishing of rational coefficients is checked with sp.cancel() == 0.
  * constants at x -> 0+ are evaluated with the proved special values
    Li2(0)=0, Li2(1)=pi^2/6, Li2(-1)=-pi^2/12, Li2(1/2)=pi^2/12-log(2)^2/2,
    Li3(0)=0, Li3(1)=zeta(3), Li3(-1)=-3 zeta(3)/4,
    Li3(1/2)=7 zeta(3)/8 - pi^2 log2 /12 + log2^3/6,
    substituted BEFORE comparison; comparison is coefficient-by-coefficient
    in the formal basis {zeta(3), pi, log(2)} (sufficient: the formal
    combination being 0 forces the real number to be 0).
"""
import json, re
import sympy as sp

x = sp.symbols('x', positive=True)
LX = sp.symbols('LX')          # formal symbol for log(x) inside series work
z3 = sp.zeta(3)
P2 = sp.pi**2
L2 = sp.log(2)
logx = sp.log(x)
l1m = sp.log(1 - x)
l1p = sp.log(1 + x)

# ---- the 8 standard arguments -------------------------------------------------
ARGS = {
    'x': x,
    '-x': -x,
    '1-x': 1 - x,
    '1/(1+x)': 1/(1 + x),
    '(1+x)/2': (1 + x)/2,
    '(1-x)/(1+x)': (1 - x)/(1 + x),
    '-(1-x)/(1+x)': -(1 - x)/(1 + x),
    '2x/(1+x)': 2*x/(1 + x),
}
ARG_NAMES = list(ARGS.keys())

# -log(1 - a(x)) expanded by hand into {log x, log(1-x), log(1+x), log 2}.
# validity: every factor below is positive for 0 < x < 1, so log(u*v) = log u + log v
# and log(u/v) = log u - log v apply literally.
NEG_LOG_1M = {
    'x':              -l1m,                       # -log(1-x)
    '-x':             -l1p,                       # -log(1+x)
    '1-x':            -logx,                      # -log(x)
    '1/(1+x)':        -(logx - l1p),              # -log(x/(1+x))
    '(1+x)/2':        -(l1m - L2),                # -log((1-x)/2)
    '(1-x)/(1+x)':    -(L2 + logx - l1p),         # -log(2x/(1+x))
    '-(1-x)/(1+x)':   -(L2 - l1p),                # -log(2/(1+x))
    '2x/(1+x)':       -(l1m - l1p),               # -log((1-x)/(1+x))
}

def arg_key(z):
    """match a sympy rational expression to one of the 8 argument names, or None."""
    zc = sp.cancel(z)
    for name in ARG_NAMES:
        if sp.cancel(zc - ARGS[name]) == 0:
            return name
    return None

def replace_polylog1(expr):
    """polylog(1, a) -> -log(1-a) via the hand-expanded table (exact on (0,1))."""
    out = expr
    for pl in sp.preorder_traversal(expr):
        if isinstance(pl, sp.polylog) and pl.args[0] == 1:
            name = arg_key(pl.args[1])
            assert name is not None, f"polylog(1, .) with unknown argument {pl.args[1]}"
            out = out.subs(pl, NEG_LOG_1M[name])
    return out

def dpoly(expr):
    """x-derivative with polylog(1,.) immediately rewritten in the log basis."""
    return replace_polylog1(sp.diff(expr, x))

# ---- parser for classical_dict.json -------------------------------------------
def parse_function_expr(name):
    if name == '1':
        return sp.Integer(1)
    out = sp.Integer(1)
    for part in name.split('*'):
        if part.startswith('Li3('):
            out *= sp.polylog(3, ARGS[part[4:-1]])
        elif part.startswith('Li2('):
            out *= sp.polylog(2, ARGS[part[4:-1]])
        elif part == 'lx':
            out *= logx
        elif part == 'l1m':
            out *= l1m
        elif part == 'l1p':
            out *= l1p
        elif part == 'pi2':
            out *= P2
        elif part == 'log2':
            out *= L2
        elif part == 'log2^2':
            out *= L2**2
        elif part == 'log2^3':
            out *= L2**3
        elif part == 'z3':
            out *= z3
        else:
            raise ValueError(f"unknown factor {part!r} in {name!r}")
    return out

def build_classical_expr(entry):
    """entry = list of [function_name, const_name, coeff_string]."""
    expr = sp.Integer(0)
    for func_name, const_name, coeff in entry:
        expr += sp.Rational(coeff) * parse_function_expr(const_name) * parse_function_expr(func_name)
    return sp.expand(expr)

def load_classical_dict(path):
    with open(path) as f:
        raw = json.load(f)
    return {w: build_classical_expr(entry) for w, entry in raw.items()}

# ---- f_a kernels ----------------------------------------------------------------
def f_a(a):
    return {0: 1/x, 1: 1/(1 - x), -1: 1/(1 + x)}[a]

# ---- special values -----------------------------------------------------------
LI2_VAL = {   # Li2 at the limit points of the 8 arguments at x = 0
    0:  sp.Integer(0),
    1:  P2/6,
    -1: -P2/12,
    sp.Rational(1, 2): P2/12 - L2**2/2,
}
LI3_VAL = {
    0:  sp.Integer(0),
    1:  z3,
    -1: -sp.Rational(3, 4)*z3,
    sp.Rational(1, 2): sp.Rational(7, 8)*z3 - P2*L2/12 + L2**3/6,
}

def _series_rat(expr, n):
    """exact Taylor/Laurent series of a rational function of x about 0, no O()."""
    s = sp.series(expr, x, 0, n + 1)
    return sp.expand(s.removeO())

def _dLi2_at(u0, k):
    """d^k/du^k Li2(u) at u0 in (-1,1)\{0}, exact (uses Li2' = -log(1-u)/u)."""
    u = sp.symbols('u')
    return sp.diff(sp.polylog(2, u), u, k).subs(u, u0)

def li2_expansion(name, n=4):
    """Li2(a(x)) about x=0 through x^n, LX = log x; valid mod O(x^{n+1} log x)."""
    a = ARGS[name]
    a0 = sp.limit(a, x, 0)
    if a0 == 0:
        s = _series_rat(a, n)                      # a = a1 x + a2 x^2 + ...
        return sp.expand(sum(s**k/sp.Integer(k)**2 for k in range(1, n + 1)))
    if a0 in (sp.Rational(1, 2), -1):
        d = _series_rat(a - a0, n)
        return sp.expand(LI2_VAL[a0] + sum(_dLi2_at(a0, k)*d**k/sp.factorial(k)
                                           for k in range(1, n + 1)))
    assert a0 == 1, name
    y = _series_rat(1 - a, n)                      # y = y1 x + y2 x^2 + ..., y1 > 0
    y1 = sp.Rational(sp.expand(y).coeff(x, 1))
    # log y = LX + log y1 + log(1 + (y/(y1 x) - 1)); expand the last log
    v = _series_rat(y/(y1*x) - 1, n - 1)           # v = v1 x + v2 x^2 + ...
    logy = LX + sp.log(y1) + sum((-1)**(k + 1)*v**k/sp.Integer(k) for k in range(1, n))
    # Li2(1-y) = pi^2/6 + sum_{k>=0} y^{k+1} (log y/(k+1) - 1/(k+1)^2)
    return sp.expand(P2/6 + sum(y**(k + 1)*(logy/sp.Integer(k + 1) - 1/sp.Integer(k + 1)**2)
                                for k in range(n)))

def li3_expansion(name, n=2):
    """Li3(a(x)) about x=0 through x^n (enough: Li3 never has a pole coefficient)."""
    a = ARGS[name]
    a0 = sp.limit(a, x, 0)
    if a0 == 0:
        s = _series_rat(a, n)
        return sp.expand(sum(s**k/sp.Integer(k)**3 for k in range(1, n + 1)))
    if a0 in (sp.Rational(1, 2), -1):
        u = sp.symbols('u')
        d = _series_rat(a - a0, n)
        return sp.expand(LI3_VAL[a0] + sum(sp.diff(sp.polylog(3, u), u, k).subs(u, a0)*d**k/sp.factorial(k)
                                           for k in range(1, n + 1)))
    assert a0 == 1, name
    # Li3(1-y) = zeta(3) - pi^2/6 y + O(y^2 log y); derived from Li3' = Li2/u.
    y = _series_rat(1 - a, 1)
    return sp.expand(z3 - P2/6*y)

def taylor0(expr, n=4):
    """Expand expr about x = 0+ through x^n (Laurent), with LX for log x.

    Returns a plain polynomial in x and LX (no O term).  Exact mod
    O(x^{n+1} log^k x).  Every polylog/log atom is replaced by its hand-coded
    expansion BEFORE expansion, so no polylog simplification is ever needed."""
    rep = {}
    for pl in sp.preorder_traversal(expr):
        if isinstance(pl, sp.polylog):
            s, a = pl.args
            name = arg_key(a)
            assert name is not None, f"polylog with unknown argument {a}"
            rep[pl] = li2_expansion(name, n) if s == 2 else li3_expansion(name, n)
    e = expr.subs(rep, simultaneous=True)
    e = e.subs({logx: LX,
                l1m: _series_rat(sp.log(1 - x), n),
                l1p: _series_rat(sp.log(1 + x), n)})
    e = sp.expand(e)
    s = sp.series(e, x, 0, n + 1)          # exact: e is now rational in x
    return sp.expand(s.removeO())

def const_term(expr, n=4):
    """coefficient of x^0 LX^0 in taylor0; also returns the residual terms."""
    s = taylor0(expr, n)
    c0 = sp.expand(s).coeff(x, 0).coeff(LX, 0)
    resid = sp.expand(s - c0)
    return sp.expand(c0), resid

def is_zero_const(c0):
    """c0 is a polynomial in zeta(3), pi, log(2) with rational coefficients;
    check it is the zero POLYNOMIAL (sufficient for the real number to be 0)."""
    p = sp.Poly(sp.expand(c0), sp.zeta(3), sp.pi, sp.log(2),
                domain=sp.QQ, extension=None)
    return all(sp.Rational(c) == 0 for c in p.coeffs()), p

# ---- log-basis coefficient extraction -------------------------------------------
LOG_ATOMS = [logx, l1m, l1p, L2]

def check_weight1_zero(expr):
    """expr must be sum of rational(x) * single log atom + rational(x).
    Returns (ok, dict atom -> cancelled coefficient)."""
    e = sp.expand(expr)
    coeffs = {}
    ok = True
    rem = e
    for atom in LOG_ATOMS:
        c = sp.cancel(e.coeff(atom))
        coeffs[str(atom)] = c
        if c != 0:
            ok = False
        rem = rem.subs(atom, 0)
    rem = sp.cancel(rem)
    coeffs['1'] = rem
    if rem != 0:
        ok = False
    return ok, coeffs
'''
with open("cor3_common.py", "w") as f:
    f.write(common)
print("written", len(common), "chars")
