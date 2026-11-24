"""Q1S: exact arithmetic in F = Q(zeta8), S-unit factorization (S = primes above 2,3).
x = (c0,c1,c2,c3) = c0 + c1*z + c2*z^2 + c3*z^3, z^4 = -1, Fractions.
V-basis: a=1+sqrt2 (unit), p2=1-zeta (N=2), q1=sqrt2+i (N=9), q2=sqrt2-i (N=9).
Facts: q1 = i*(1 - i*sqrt2), q2 = -i*(1 + i*sqrt2); sigma: z->z^3 fixes sqrt(-2)=z+z^3;
N_{F/Q(sqrt-2)}(x) = x*sigma(x);  v_{q1}(x) = v_{(1-sqm2)}(N x), v_{q2}(x) = v_{(1+sqm2)}(N x).
Torsion <zeta> order 8. sqrt2 = z - z^3, i = z^2.
"""
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp

ZERO=(Fr(0),)*4; ONE=(Fr(1),Fr(0),Fr(0),Fr(0)); ZETA=(Fr(0),Fr(1),Fr(0),Fr(0))
def fadd(x,y): return tuple(a+b for a,b in zip(x,y))
def fneg(x): return tuple(-a for a in x)
def fsub(x,y): return tuple(a-b for a,b in zip(x,y))
def fmul(x,y):
    r=[Fr(0)]*7
    for i in range(4):
        for j in range(4): r[i+j]+=x[i]*y[j]
    return tuple(r[k]-(r[k+4] if k+4<7 else Fr(0)) for k in range(4))
def fpow(x,k):
    if k<0: return fpow(finv(x),-k)
    r=ONE
    for _ in range(k): r=fmul(r,x)
    return r
def finv(x):
    n=fnorm(x); assert n!=0
    others=ONE
    for j in (3,5,7): others=fmul(others,fconj_z(x,j))
    return fscale(others, Fr(1)/n)
def fconj_z(x,j):
    zj=fpow(ZETA,j%8); r=ZERO
    for k in range(4): r=fadd(r,fscale(fpow(zj,k),x[k]))
    return r
def fconj(x): return fconj_z(x,7)
def fsigma(x): return fconj_z(x,3)   # fixes sqrt(-2)
def fscale(x,q):
    q=Fr(q); return tuple(q*a for a in x)
def fnorm(x):
    n=ONE
    for j in (1,3,5,7): n=fmul(n,fconj_z(x,j))
    assert n[1]==0 and n[2]==0 and n[3]==0
    return n[0]
def feq(x,y): return fsub(x,y)==ZERO
def fshow(x): return "(%s)+(%s)z+(%s)z^2+(%s)z^3"%tuple(str(c) for c in x)
def femb(x,j):
    w=mp.exp(mp.pi*1j*j/4)
    return sum(mp.mpf(c.numerator)/c.denominator*w**k for k,c in enumerate(x))

z=ZETA; i_=fpow(ZETA,2); r2=fsub(ZETA,fpow(ZETA,3))
a_=fadd(ONE,r2); p2_=fsub(ONE,ZETA); q1_=fadd(r2,i_); q2_=fsub(r2,i_)
BASISEL={'a':a_,'p2':p2_,'q1':q1_,'q2':q2_}; BASIS=('a','p2','q1','q2')
sqm2=fadd(ZETA,fpow(ZETA,3))       # sqrt(-2) = z + z^3

def _vpq(R,S,sign):
    """valuation of R+S*sqm2 (integers) at prime (1+sign*sqm2) of Z[sqrt-2]."""
    e=0
    while True:
        # (R + S m)/(1+sign m) = (R + 2 sign S)/3 + (S - sign R)/3 m
        nr = R + 2*sign*S; ns = S - sign*R
        if nr%3==0 and ns%3==0: R,S = nr//3, ns//3; e+=1
        else: return e

def factor_sunit(x):
    """-> (k, dict) with x = zeta^k * a^ea * p2^e2 * q1^e1 * q2^e2p, or None if not an S-unit.
    Vector for the symbol: dict over BASIS (torsion dropped)."""
    if feq(x,ZERO): return None
    N=fnorm(x); num,den=abs(N.numerator),abs(N.denominator)
    def vp(n,p):
        e=0
        while n%p==0: n//=p; e+=1
        return e,n
    e2n,num1=vp(num,2); e3n,num2=vp(num1,3)
    e2d,den1=vp(den,2); e3d,den2=vp(den1,3)
    if num2!=1 or den2!=1: return None
    e2=e2n-e2d; B=e3n-e3d
    # split q1/q2 via relative norm
    w=fmul(x,fsigma(x))
    assert w[2]==0 and w[1]==w[3]
    r,s=w[0],w[1]
    den=r.denominator*s.denominator//gcd(r.denominator,s.denominator)
    R=r.numerator*(den//r.denominator); S=s.numerator*(den//s.denominator)
    v3=0; d=den
    while d%3==0: d//=3; v3+=1
    v_q  = _vpq(R,S,+1) - v3     # v_{(1+sqm2)}(w);  w = N_{F/k}(x)
    v_qb = _vpq(R,S,-1) - v3     # v_{(1-sqm2)}(w)
    # primes above 3 are inert in F/k (residue degree 2): v_q(w) = 2 v_{q2}(x), v_qb(w) = 2 v_{q1}(x)
    assert v_q%2==0 and v_qb%2==0
    e1p, e1 = v_q//2, v_qb//2
    assert 2*(e1+e1p)==B, (x,e1,e1p,B)
    # unit exponent e_a via embedding 1: |sigma1(u)| = a^ea
    u=fmul(x, fmul(fpow(p2_,-e2), fmul(fpow(q1_,-e1), fpow(q2_,-e1p))))
    with mp.workdps(50):
        au=abs(femb(u,1)); la=mp.log(1+mp.sqrt(2))
        ea=int(mp.nint(mp.log(au)/la))
        # torsion k from argument
        k=int(mp.nint(4*mp.arg(femb(u,1))/mp.pi))%8
    # exact verify
    rec=fmul(fpow(ZETA,k), fmul(fpow(a_,ea), fmul(fpow(p2_,e2), fmul(fpow(q1_,e1), fpow(q2_,e1p)))))
    assert feq(rec,x), ("factorization failed", x, (k,ea,e2,e1,e1p))
    vec={}
    if ea: vec['a']=Fr(ea)
    if e2: vec['p2']=Fr(e2)
    if e1: vec['q1']=Fr(e1)
    if e1p: vec['q2']=Fr(e1p)
    return (k, vec)

if __name__=="__main__":
    # self-tests
    assert fnorm(p2_)==2, fnorm(p2_)
    assert fnorm(q1_)==9 and fnorm(q2_)==9
    assert fnorm(a_)==1
    # q1 = i(1 - i sqrt2):
    assert feq(q1_, fmul(i_, fsub(ONE, fmul(i_,r2))))
    assert feq(q2_, fmul(fneg(i_), fadd(ONE, fmul(i_,r2))))
    for nm,x in [('a',a_),('p2',p2_),('q1',q1_),('q2',q2_),('sqrt2',r2),('i',i_),('1+i',fadd(ONE,i_))]:
        print(nm, factor_sunit(x))
