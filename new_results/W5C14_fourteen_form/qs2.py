"""Exact arithmetic in Q(sqrt2): elements (p,q) = p+q*sqrt2, Fractions."""
from fractions import Fraction as F
class Q:
    __slots__=('p','q')
    def __init__(s,p=0,q=0): s.p=F(p); s.q=F(q)
    def __add__(s,o): o=_c(o); return Q(s.p+o.p,s.q+o.q)
    __radd__=__add__
    def __neg__(s): return Q(-s.p,-s.q)
    def __sub__(s,o): return s+(-_c(o))
    def __rsub__(s,o): return _c(o)-s
    def __mul__(s,o): o=_c(o); return Q(s.p*o.p+2*s.q*o.q, s.p*o.q+s.q*o.p)
    __rmul__=__mul__
    def conj(s): return Q(s.p,-s.q)
    def norm(s): return s.p*s.p-2*s.q*s.q
    def inv(s): n=s.norm(); return Q(s.p/n,-s.q/n)
    def __truediv__(s,o): return s*_c(o).inv()
    def __rtruediv__(s,o): return _c(o)*s.inv()
    def __eq__(s,o): o=_c(o); return s.p==o.p and s.q==o.q
    def __hash__(s): return hash((s.p,s.q))
    def iszero(s): return s.p==0 and s.q==0
    def __float__(s): return float(s.p)+float(s.q)*2**0.5
    def __repr__(s): return f'({s.p}{"+" if s.q>=0 else ""}{s.q}*r2)'
def _c(o): return o if isinstance(o,Q) else Q(o,0)
R2=Q(0,1); A=Q(1,1); B=Q(-1,1)
