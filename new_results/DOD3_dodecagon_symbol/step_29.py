import os
code = '''"""DOD3: exact symbolic R(i;j,k) = T3a - (SwG-SwD)/x_i for all 12 triples.
Adapted from DOD2 step_16 with EXACT logfac via S-unit factorization (no PSLQ).
Logs: l2=log2, l3=log3, lu=log(2-sqrt3). For S-unit theta = +-p2^v2 p3^v3 eps^k:
log|theta| = (v2/2) l2 + (v3/2) l3 + (-v2/2 - k) lu.   [log p2=(l2-lu)/2, log p3=l3/2, log eps=-lu]
"""
import sympy as sp
from sympy import Rational, sqrt, symbols, Function, expand, Add
from fractions import Fraction as Fr
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dod3_symbol import *

r3 = sqrt(3)
l2,l3,lu = symbols('l2 l3 lu')
Li2A = Function('Li2A'); Li3A = Function('Li3A')

def sp_of(a):
    m,n,d = a
    return Rational(m,d) + Rational(n,d)*r3
def pair_of(e):
    """sympy e in Q(sqrt3) -> pair (m,n,d)"""
    e = sp.expand(e)
    p = sp.simplify((e + e.subs(r3,-r3))/2)
    q = sp.simplify((e - e.subs(r3,-r3))/(2*r3))
    p = Rational(p); q = Rational(q)
    from math import lcm
    d = lcm(p.denominator, q.denominator)
    return P(p*d, q*d, d)

Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
IDX = [1,2,5,6]
Xsym = {i: sp_of(Xs[i]) for i in IDX}

def logfac_pair(th):
    """th: pair != 0 -> symbolic log|th|"""
    if peq(th, P(0,0)): raise ValueError('log 0')
    sgn = 1
    if pfloat(th) < 0: th = pneg(th)
    if peq(th, ONE): return sp.Integer(0)
    _,v2,v3,k = factor_sunit(th)
    return Fr(v2,2)*l2 + Fr(v3,2)*l3 + (-Fr(v2,2)-k)*lu
def logfac_sym(the):
    return logfac_pair(pair_of(the))

def L1abs(kap, u):
    """kap, u pairs. integral int log|1-kap t|/t type piece."""
    ku = pmul(kap,u)
    if pfloat(psub(ku,ONE)) <= 0:
        return -Li2A(sp_of(ku))
    else:
        return logfac_pair(kap)*logfac_pair(u) + logfac_pair(u)**2/2 + Li2A(sp_of(pdiv(ONE,ku)))
def I2abs(kap, u):
    ku = pmul(kap,u)
    if pfloat(psub(ku,ONE)) <= 0:
        return -logfac_pair(u)*Li2A(sp_of(ku)) + Li3A(sp_of(ku))
    else:
        return (logfac_pair(kap)*logfac_pair(u)**2/2 + logfac_pair(u)**3/3
                + logfac_pair(u)*Li2A(sp_of(pdiv(ONE,ku))) + Li3A(sp_of(pdiv(ONE,ku))))
def Mcf(c, xi):
    """c, xi pairs"""
    kap = pdiv(c, psub(c,ONE))
    lo, hi = psub(ONE,xi), ONE
    return ( logfac_pair(pneg(psub(ONE,c)) if pfloat(psub(ONE,c))<0 else psub(ONE,c))*(logfac_pair(hi)**2-logfac_pair(lo)**2)/2
             + I2abs(kap,hi) - I2abs(kap,lo)
             - (logfac_pair(hi)**3-logfac_pair(lo)**3)/3 )
def QnonC(a,b,xi):
    ka = pdiv(a, psub(a,ONE)); kb = pdiv(b, psub(b,ONE))
    A = logfac_pair(psub(ONE,a) if pfloat(psub(ONE,a))>0 else pneg(psub(ONE,a)))
    B = logfac_pair(psub(ONE,b) if pfloat(psub(ONE,b))>0 else pneg(psub(ONE,b)))
    lo, hi = psub(ONE,xi), ONE
    def F(u):
        return ( A*B*logfac_pair(u) + A*L1abs(kb,u) - A*logfac_pair(u)**2/2
               + B*L1abs(ka,u) - B*logfac_pair(u)**2/2
               - I2abs(ka,u) - I2abs(kb,u) + logfac_pair(u)**3/3 )
    return F(hi)-F(lo)
def Swrem(ka,kb):
    A = logfac_pair(psub(ONE,ka) if pfloat(psub(ONE,ka))>0 else pneg(psub(ONE,ka)))
    return -A*(Li2A(sp_of(kb))-Li2A(sp_of(pdiv(kb,ka)))) + Li3A(sp_of(kb))

def build_R(i,j,k):
    xi = Xs[i]
    a = pdiv(Xs[j],xi); b = pdiv(Xs[k],xi)
    ka = pdiv(a, psub(a,ONE)); kb = pdiv(b, psub(b,ONE))
    lo = psub(ONE,xi)
    da, db = pmul(ka,lo), pmul(kb,lo)
    bnd = -logfac_pair(lo)**3/3
    R = (bnd + Mcf(a,xi) + Mcf(b,xi) + QnonC(a,b,xi) + Swrem(ka,kb) - Swrem(da,db))
    # divide by xi (pair): multiply sympy expr by sp_of(1/xi)
    R = sp.expand(R * sp_of(pdiv(ONE,xi)))
    return R

TRIPLES = [(1,2,5),(1,2,6),(1,5,6),(2,1,5),(2,1,6),(2,5,6),(5,1,2),(5,1,6),(5,2,6),(6,1,2),(6,1,5),(6,2,5)]
'''
with open('dod3_R.py','w') as f:
    f.write(code)
print("wrote dod3_R.py")
# quick self-test on (2,5,6)
import sys; sys.path.insert(0, os.getcwd())
from dod3_R import *
R = build_R(2,5,6)
print("R(2,5,6) built, #terms:", len(Add.make_args(R)))
atoms = R.atoms(Li2A) | R.atoms(Li3A)
print("polylog atoms:", sorted(str(s) for s in atoms))
