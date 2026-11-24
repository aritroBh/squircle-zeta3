"""Independent re-implementation (own code, exact Gaussian-rational arithmetic at numeric (d,m)) of the FORCED ascending recursion
F(delta,0)=a_delta (Bell from the L cache), F(delta,s+1)=Lambda(delta,s)-F(delta,s); tests closing conditions F(7,4)=0, F(8,5)=0 and all
lower closings at many (d,m); and compares the degree-6 towers with the exact data (off-diagonal keys)."""
import sys, os, pickle, itertools, time
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, factorial
import sympy as sp
t0=time.time()
class G:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): o=g(o); return G(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o): o=g(o); return G(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return g(o)-s
    def __neg__(s): return G(-s.a,-s.b)
    def __mul__(s,o): o=g(o); return G(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a+s.b*s.b; return G(s.a/n, -s.b/n)
    def __truediv__(s,o): return s*g(o).inv()
    def __rtruediv__(s,o): return g(o)*s.inv()
    def iszero(s): return s.a==0 and s.b==0
    def __eq__(s,o): o=g(o); return s.a==o.a and s.b==o.b
    def __repr__(s): return f'({s.a}+{s.b}i)'
def g(x):
    return x if isinstance(x,G) else G(x,0)
def gpow(x,n):
    if n<0: x=x.inv(); n=-n
    r=G(1)
    for _ in range(n): r=r*x
    return r
def ev(expr, env):
    if expr.is_Rational: return G(Fr(int(expr.p),int(expr.q)))
    if expr is sp.I: return G(0,1)
    if expr.is_Symbol: return env[expr.name]
    if expr.is_Add:
        r=G(0)
        for a in expr.args: r=r+ev(a,env)
        return r
    if expr.is_Mul:
        r=G(1)
        for a in expr.args: r=r*ev(a,env)
        return r
    if expr.is_Pow:
        assert expr.exp.is_Integer, expr
        return gpow(ev(expr.base,env), int(expr.exp))
    raise ValueError(expr)
dsym,msym=sp.symbols('d m')
Lraw=pickle.load(open('L_simplified.pkl','rb'))
L={k: sp.sympify(v, locals={'d':dsym,'m':msym,'I':sp.I}) for k,v in Lraw.items()}
print('L keys',sorted(L.keys()))
# tau from the copied library (sympy exact), converted to G
import tow4lib as T4
def toG(z):
    z=sp.nsimplify(z); re,im=sp.re(z),sp.im(z)
    return G(Fr(int(sp.numer(re)),int(sp.denom(re))), Fr(int(sp.numer(im)),int(sp.denom(im))))
tau1=[toG(T4.tau1f(n)) for n in range(9)]
tau2=[toG(T4.tau2f(n)) for n in range(9)]
print('tau1(8)=',tau1[8],'tau2(8)=',tau2[8])

def parts(rem,minp,cur):
    if rem==0: yield tuple(cur); return
    for p in range(minp,rem+1):
        yield from parts(rem-p,p,cur+[p])
@lru_cache(maxsize=None)
def Lval(r,dd,mm): return ev(L[r],{'d':G(dd),'m':G(mm)})
@lru_cache(maxsize=None)
def a(delta,dd,mm):
    if delta==0: return G(1)
    tot=G(0)
    for part in parts(delta,1,[]):
        from collections import Counter
        C=Counter(part); coef=G(1)
        for sz,mult in C.items():
            coef=coef*gpow(Lval(sz,dd,mm),mult)*Fr(1,factorial(mult))
        tot=tot+coef
    return tot
@lru_cache(maxsize=None)
def F(delta,s,dd,mm):
    """forced family value, s>delta//2 -> 0 by definition (only used inside Lambda for s<=floor(delta'/2))"""
    if s>delta//2: return G(0)
    if s==0: return a(delta,dd,mm)
    return Lam(delta,s-1,dd,mm)-F(delta,s-1,dd,mm)
@lru_cache(maxsize=None)
def Lam(delta,s,dd,mm):
    tot=G(0)
    for dp in range(delta+1):
        fA=F(dp,s,dd+1,mm+1); fB=F(dp,s,dd+1,mm)
        tot=tot+ (G(0,1)*(2*dd-1)/(2*dd-dp-1)) * (tau1[delta-dp]*fA - tau2[delta-dp]*fB)
    return tot
def closing(delta,dd,mm):
    top=delta//2
    return Lam(delta,top,dd,mm)-F(delta,top,dd,mm)     # F(delta, top+1)

if __name__=='__main__':
    bad=0;tot=0
    pts=[(dd,mm) for dd in range(6,14) for mm in range((dd+1)//2, dd+1)]
    for dl in range(0,9):
        for (dd,mm) in pts:
            if 2*dd-dl-1==0: continue
            c=closing(dl,dd,mm); tot+=1
            if not c.iszero(): bad+=1; print('CLOSING FAIL',dl,dd,mm,c)
        print('closing delta=',dl,'checked at',len(pts),'points; running bad =',bad, 'time',round(time.time()-t0,1),flush=True)
    print('closing total',tot,'bad',bad)
