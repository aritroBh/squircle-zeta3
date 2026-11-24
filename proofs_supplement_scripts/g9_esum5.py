# Copied from research/08_theory/K1f_hexagon_order5/esum5.py (our notes) with the working directory taken from the environment variable G9WORK.
# Only paths and the module name of the e-sum engine were changed.
import os

# esum5.py: weight-5 e-sum engine. Reduces S = sum_{k>=K} x^k g(k) to the HPL atom basis,
# where g(k) = sum of products of harmonic numbers Hf(j,k+d)^{e} and rational functions of k
# with poles at integers (poles at +d must be < K). Total weight must be <= 5.
# Atoms from atom_table.pkl (verified in step_26). x kept symbolic.
import sympy as sp
from functools import lru_cache

x = sp.Symbol('x')
H1,H2,H3,H4 = sp.symbols('H1 H2 H3 H4')     # harmonic symbols at the current base variable
HSYMS = {1:H1, 2:H2, 3:H3, 4:H4}
Hf = sp.Function('Hf')                       # Hf(j, m) = H_m^{(j)}

_atom_tab = None
def load_atom_table(WD):
    global _atom_tab
    raw = __import__('pickle').load(open(f'{WD}/atom_table.pkl','rb'))
    _atom_tab = {k: sp.sympify(v) for k, v in raw.items()}

def atom_expr(a, ms):
    return _atom_tab[(a, tuple(sorted(ms)))]

def shift_h(g, k):
    """align all Hf(j, m) args to base variable k."""
    def rep(f):
        j = int(f.args[0]); m = f.args[1]
        d = sp.expand(m - k)
        if not d.is_integer: raise ValueError(f"non-integer shift {f}")
        d = int(d)
        if d == 0: return HSYMS[j]
        if d > 0: return HSYMS[j] + sum(1/(k+i)**j for i in range(1, d+1))
        return HSYMS[j] - sum(1/(k-i)**j for i in range(0, -d))
    return sp.expand(g.replace(lambda f: isinstance(f, sp.Function) and f.func==Hf, rep))

def splitH(term, k):
    """term -> (e_tuple (4,), R(k)) with e_j the exponents of HSYMS."""
    e = [0,0,0,0]
    R = term
    for j in range(1,5):
        hj = HSYMS[j]
        if term.has(hj):
            pj = sp.degree(term, hj)
            e[j-1] = pj
            R = sp.expand(R / hj**pj)
    return tuple(e), R

def pieces(R, k):
    """partial fractions of R in k -> (c, kind, d, a): kind 'z' pole 0, 'p' pole -d, 'm' pole +d."""
    if R == 0: return
    for piece in sp.Add.make_args(sp.apart(R, k)):
        fs = sp.Mul.make_args(piece) if isinstance(piece, sp.Mul) else [piece]
        c = sp.Integer(1); den = sp.Integer(1)
        for f in fs:
            if f.has(k): den *= f
            else: c *= f
        if den == 1:
            raise ValueError(f"polynomial part {piece} in {R}")
        base_, exp_ = den.as_base_exp()
        a_ = int(-exp_); assert a_ >= 1
        if base_ == k:
            yield c, 'z', 0, a_
        else:
            d = int(sp.expand(base_ - k))
            if d >= 1: yield c, 'p', d, a_
            else: yield c, 'm', -d, a_

def Ssum(g, K, k):
    """sum_{k>=K} x^k g(k); g = poly in HSYMS with rational coeffs having poles at -s (s>=0) or +d (1<=d<K)."""
    tot = sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(g)):
        e, R = splitH(term, k)
        for c, kind, d, a in pieces(R, k):
            if kind in ('z','p'):
                v = base(e, d, a, k)
                for n0 in range(1, K):
                    v -= x**n0 * sp.prod(sp.harmonic(n0, j)**ej for j, ej in enumerate(e, start=1) if ej) / (n0+d)**a
                tot += c*v
            else:
                assert d < K, (R, K)
                # substitute k = m + d: Hf aligned then recurse
                g2 = c*sp.prod(Hf(j, k+d)**ej for j, ej in enumerate(e, start=1) if ej) / k**a
                g2 = shift_h(g2, k)   # now in terms of HSYMS at k... careful: this replaces k as var
                tot += sp.expand(x**d * Ssum(g2, K-d, k))
    return sp.expand(tot)

@lru_cache(None)
def base_cached(e, s, a):
    return _base(e, s, a)

def base(e, s, a, k=None):
    """sum_{m>=1} x^m prod (H_m^{(j)})^e_j /(m+s)^a, s>=0."""
    if s == 0:
        ms = []
        for j, ej in enumerate(e, start=1): ms += [j]*ej
        return atom_expr(a, ms)
    m = sp.Symbol('m')
    g = sp.prod((Hf(j, m) - sum(1/(m-i)**j for i in range(0, s)))**ej
                for j, ej in enumerate(e, start=1) if ej) / m**a
    g = shift_h(g, m)
    return sp.expand(x**(-s) * Ssum(g, s+1, m))

def esum(g, K, k):
    """user entry: sum_{k>=K} x^k g(k), g may contain Hf(j, m) with m affine in k."""
    return Ssum(shift_h(g, k), K, k)
