"""DOD3: symbol machinery over F = Q(sqrt3), S = {p2=1+sqrt3, p3=sqrt3}, units {+-1, eps=2+sqrt3}.
V = O_{F,S}^x / tors (x) Q = Q^3 with basis ('2','3','e').
Symbols in V^{otimes3} = Q^27, dict representation.
Exact Q(sqrt3) elements as pairs (m,n,d) = (m+n sqrt3)/d.
"""
from fractions import Fraction as Fr
from math import gcd
import math, itertools
import sympy as sp

# ---------- exact Q(sqrt3) ----------
def P(m, n, d=1):
    m, n, d = int(m), int(n), int(d)
    if d < 0: m, n, d = -m, -n, -d
    g = gcd(gcd(abs(m), abs(n)), d)
    return (m//g, n//g, d//g)
def padd(a,b):
    m1,n1,d1=a; m2,n2,d2=b
    return P(m1*d2+m2*d1, n1*d2+n2*d1, d1*d2)
def pneg(a): m,n,d=a; return P(-m,-n,d)
def psub(a,b): return padd(a,pneg(b))
def pmul(a,b):
    m1,n1,d1=a; m2,n2,d2=b
    return P(m1*m2+3*n1*n2, m1*n2+n1*m2, d1*d2)
def pdiv(a,b):
    m2,n2,d2=b
    den = m2*m2-3*n2*n2
    return pmul(a, P(m2*d2,-n2*d2,den))
def ppow(a,k):
    if k<0: return ppow(pdiv(P(1,0),a),-k)
    r=P(1,0)
    for _ in range(k): r=pmul(r,a)
    return r
def pscale(a,q):
    m,n,d=a
    return P(m*q.numerator, n*q.numerator, d*q.denominator)
def pfloat(a):
    m,n,d=a; return (m+n*(3.0**0.5))/d
def peq(a,b): return a==b
def psplit(a):
    m,n,d=a; return Fr(m,d), Fr(n,d)
def pconj(a):
    m,n,d=a; return P(m,-n,d)

ONE=P(1,0); R3=P(0,1)
P2=P(1,1); P3=P(0,1); EPS=P(2,1)
LN_EPS = math.log(2+3**0.5)
def pnorm(a):
    m,n,d=a; return Fr(m*m-3*n*n, d*d)
def vp23(x):
    e2=e3=0
    while x%2==0: x//=2; e2+=1
    while x%3==0: x//=3; e3+=1
    return e2,e3,x
def factor_sunit(a):
    """a = sign * P2^v2 * P3^v3 * EPS^k  or ('NOT_SUNIT', info)"""
    if peq(a,P(0,0)): return None
    N=pnorm(a); num,den=abs(N.numerator),N.denominator
    e2n,e3n,rn=vp23(num); e2d,e3d,rd=vp23(den)
    if rn!=1 or rd!=1: return ('NOT_SUNIT',('residual',rn,rd))
    v2=e2n-e2d; v3=e3n-e3d
    u=pmul(a,pmul(ppow(P2,-v2),ppow(P3,-v3)))
    fu=abs(pfloat(u))
    k=int(round(math.log(fu)/LN_EPS)) if fu>0 else 0
    for kk in (k-1,k,k+1):
        if peq(u,ppow(EPS,kk)): return (1,v2,v3,kk)
        if peq(u,pneg(ppow(EPS,kk))): return (-1,v2,v3,kk)
    return ('NOT_SUNIT',('unit',u))

BASIS = ('2','3','e')   # p2, p3, eps
def vecfun(a):
    """S-unit a -> vec dict over BASIS (torsion -1 dropped)."""
    f = factor_sunit(a)
    assert f is not None and f[0] != 'NOT_SUNIT', (a,f)
    _,v2,v3,k = f
    out = {}
    if v2: out['2']=Fr(v2)
    if v3: out['3']=Fr(v3)
    if k:  out['e']=Fr(k)
    return out

# ---------- vectors & symbols (over BASIS letters) ----------
def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def vmul(v,a):
    a=Fr(a); return {k:c*a for k,c in v.items() if c*a!=0}
def tensor3(u,v,w,coeff=Fr(1)):
    S={}
    for a,ca in u.items():
        for b,cb in v.items():
            for c,cc in w.items():
                S[(a,b,c)]=S.get((a,b,c),Fr(0))+coeff*ca*cb*cc
    return {k:x for k,x in S.items() if x!=0}
def sadd(*Ss):
    r={}
    for S in Ss:
        for k,c in S.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def smul(S,a):
    a=Fr(a); return {k:c*a for k,c in S.items() if c*a!=0}

def S_Li3(z):
    """symbol of Li_3(z): -(1-z) x z x z ; z an S-unit pair with 1-z S-unit."""
    zv = vecfun(z); omzv = vecfun(psub(ONE,z))
    return tensor3(vmul(omzv,-1), zv, zv)

def S_Li21(x, y):
    """symbol of Li_{2,1}(x,y), standard convention sum_{m<n} x^m y^n/(m^2 n)."""
    L1 = vecfun(psub(pmul(x,y),ONE))
    L2 = vecfun(pmul(x,y))
    L3 = vecfun(psub(y,ONE))
    L4 = vecfun(pdiv(pmul(x,psub(y,ONE)), psub(x,ONE)))
    L5 = vecfun(psub(x,ONE))
    L6 = vecfun(x)
    return sadd(tensor3(L1,L2,L3), tensor3(L1,L4,L6), tensor3(L3,L5,L6))

# ---------- quotient Q = V^3 / Decomp ----------
class Quotient:
    def __init__(self):
        self.letters = list(BASIS)
        self.dim3 = 27
        self.basis3 = list(itertools.product(self.letters, repeat=3))
        self.bidx = {t:i for i,t in enumerate(self.basis3)}
        D=[]
        for a in self.letters:
            for b in self.letters:
                for c in self.letters:
                    v = sp.zeros(1,27)
                    for t in [(a,b,c),(b,a,c),(b,c,a)]:
                        v[self.bidx[t]] += 1
                    D.append(v)
        self.D = D
        DM = sp.Matrix.hstack(*[d.T for d in D])
        self.rank_D = DM.rank()
    def coord(self,S):
        v = sp.zeros(1,27)
        for t,c in S.items():
            v[self.bidx[t]] += sp.Rational(c.numerator,c.denominator)
        return v
    def solve(self, target, gens):
        """target in span(gens)+Decomp? -> (bool, coeffs over gens)"""
        M = sp.Matrix.hstack(*[g.T for g in gens]+[d.T for d in self.D])
        r0 = M.rank()
        M2 = sp.Matrix.hstack(M, target.T)
        if M2.rank()!=r0: return False, None
        c = M.gauss_jordan_solve(target.T)[0]
        return True, [sp.Rational(x) for x in c[:len(gens)]]
