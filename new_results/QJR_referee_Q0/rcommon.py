# Referee QJR helpers (own code). Exact arithmetic in Q(sqrt2); exact logs of {2,3}-units.
import sympy as sp, mpmath as mp
r2=sp.sqrt(2); a=1+r2; b=r2-1
lam,L,l3,pi=sp.symbols('lambda L ell3 pi_',real=True)   # lambda=log2, L=log a, ell3=log3, pi_=pi
alpha=sp.Symbol('alpha_c',real=True)
def canon(x):
    x=sp.expand(sp.radsimp(sp.together(sp.sympify(x))))
    return x
def conj(x):
    return canon(sp.sympify(x).subs(r2,-r2)) if False else canon(canon(x).subs(sp.sqrt(2),-sp.sqrt(2)))
def logS(x):
    """exact log of a positive {2,3}-unit x of Q(sqrt2) as i*lam/2+j*l3+k*L; verified exactly."""
    x=canon(x); assert sp.N(x,50)>0, x
    N=canon(x*conj(x)); N=sp.Rational(N); s=abs(N)
    i=0
    while s.p%2==0: s/=2; i+=1
    while s.q%2==0: s*=2; i-=1
    j2=0
    while s.p%3==0: s/=3; j2+=1
    while s.q%3==0: s*=3; j2-=1
    assert s==1 and j2%2==0, ('not an S-unit',x,N)
    j=j2//2
    mp.mp.dps=60
    k=(mp.log(abs(mp.mpf(sp.N(x,60))))-mp.log(abs(mp.mpf(sp.N(conj(x),60)))))/(2*mp.log(1+mp.sqrt(2)))
    k=int(mp.nint(k))
    chk=canon(x/(r2**i*sp.Integer(3)**j*a**k))
    assert chk==1, ('log decomposition failed',x,i,j,k,chk)
    return sp.Rational(i,2)*lam+j*l3+k*L
Li2s={}
def Li2(x):
    x=canon(x)
    if x==0: return sp.Integer(0)
    key=str(x)
    if key not in Li2s: Li2s[key]=(sp.Function('Li2')(x))
    return Li2s[key]
def numval(expr,dps=110):
    mp.mp.dps=dps+10
    sub={lam:mp.log(2),L:mp.log(1+mp.sqrt(2)),l3:mp.log(3),pi:mp.pi,alpha:mp.atan(1/mp.sqrt(2))}
    def ev(e):
        e=sp.sympify(e)
        if e.is_Number: return mp.mpf(sp.N(e,dps+20)) if not e.is_Rational else mp.mpf(e.p)/e.q
        if e in sub: return sub[e]
        if isinstance(e,sp.Function) and e.func.__name__=='Li2': return mp.polylog(2,mp.mpf(sp.N(e.args[0],dps+20)))
        if e.is_Add: return mp.fsum(ev(t) for t in e.args)
        if e.is_Mul:
            p=mp.mpf(1)
            for t in e.args: p*=ev(t)
            return p
        if e.is_Pow: return ev(e.base)**ev(e.exp)
        if e==sp.sqrt(2): return mp.sqrt(2)
        raise ValueError(e)
    return ev(sp.expand(expr))
Lam=sp.Symbol('Lambda'); M=sp.Symbol('M')   # Lambda=Li2(b), M=Li2(-b^4)
