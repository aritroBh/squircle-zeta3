
import sympy as sp
from sympy import Rational, symbols, expand, apart, harmonic
from functools import lru_cache

k, x = sp.symbols('k x')
H, H2 = sp.symbols('H H2')
L1, L2, L3, AA, z2 = sp.symbols('L1 L2 L3 AA z2')
Hf = sp.Function('Hf'); H2f = sp.Function('H2f')

def atom(p, q, a):
    if p==0 and q==0:
        return {1: L1, 2: L2, 3: L3}[a]
    if (p,q)==(1,0) and a==1: return L2 + L1**2/sp.Integer(2)
    if (p,q)==(1,0) and a==2: return AA
    if (p,q)==(0,1) and a==1: return 3*L3 - 2*AA + L2*L1
    if (p,q)==(2,0) and a==1: return L3 + L2*L1 + L1**3/sp.Integer(3)
    raise ValueError(f"atom {(p,q,a)} not in table (weight>3?)")

@lru_cache(None)
def base(p, q, s, a):
    """sum_{n>=1} x^n H(n)^p H2(n)^q / (n+s)^a , s>=0."""
    if s == 0:
        return atom(p, q, a)
    mvar = sp.Symbol('m')
    Ht = H - sum(1/(mvar-j) for j in range(0, s))
    H2t = H2 - sum(1/(mvar-j)**2 for j in range(0, s))
    expr = sp.expand(Ht**p * H2t**q / mvar**a)
    return sp.expand(x**(-s) * Ssum(expr, s+1, mvar))

def _splitH(term):
    p_ = sp.degree(term, H) if term.has(H) else 0
    q_ = sp.degree(term, H2) if term.has(H2) else 0
    R = sp.simplify(term / (H**p_ * H2**q_))
    return p_, q_, sp.expand(R)

def _pieces(R, var):
    for piece in sp.Add.make_args(sp.apart(R, var)):
        fs = sp.Mul.make_args(piece) if isinstance(piece, sp.Mul) else [piece]
        c = sp.Integer(1); den = sp.Integer(1)
        for f in fs:
            if f.has(var): den *= f
            else: c *= f
        if den == 1:
            raise ValueError(f"polynomial part {piece}")
        base_, exp_ = den.as_base_exp()
        a_ = int(-exp_); assert a_ >= 1
        if base_ == var:
            yield c, 'z', 0, a_
        else:
            d = int(sp.expand(base_ - var))
            if d >= 1: yield c, 'p', d, a_
            else: yield c, 'm', -d, a_

def Ssum(expr, K, var):
    """sum_{n>=K} x^n expr(n), expr = sum of H^p H2^q * R(n), R's poles at -s (s>=0) or j in 1..K-1."""
    tot = 0
    for term in sp.Add.make_args(sp.expand(expr)):
        p_, q_, R = _splitH(term)
        for c, kind, j, a in _pieces(R, var):
            if kind in ('z', 'p'):
                v = base(p_, q_, j, a)
                for n0 in range(1, K):
                    v -= x**n0 * harmonic(n0)**p_ * harmonic(n0,2)**q_ / (n0+j)**a
                tot += c*v
            else:
                assert j < K, (R, K)
                nvar = sp.Symbol('n')
                Ht = H + sum(1/(nvar+jj) for jj in range(1, j+1))
                H2t = H2 + sum(1/(nvar+jj)**2 for jj in range(1, j+1))
                sub = Ssum(c*Ht**p_ * H2t**q_ / nvar**a, K-j, nvar)
                tot += x**j * sub
    return sp.expand(tot)

def shift_h(g):
    def repH(f):
        d = int(sp.expand(f.args[0] - k))
        if d == 0: return H
        if d > 0: return H + sum(1/(k+j) for j in range(1, d+1))
        return H - sum(1/(k-j+1) for j in range(1, -d+1))
    def repH2(f):
        d = int(sp.expand(f.args[0] - k))
        if d == 0: return H2
        if d > 0: return H2 + sum(1/(k+j)**2 for j in range(1, d+1))
        return H2 - sum(1/(k-j+1)**2 for j in range(1, -d+1))
    g = g.replace(lambda f: isinstance(f, sp.Function) and f.func==Hf, repH)
    g = g.replace(lambda f: isinstance(f, sp.Function) and f.func==H2f, repH2)
    return sp.expand(g)

def esum(g, K=1):
    g = shift_h(g)
    return Ssum(g, K, k)
