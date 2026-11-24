import sympy as sp, pickle, os
from fractions import Fraction as Fr
from math import comb, factorial
from sympy.functions.combinatorial.numbers import stirling
# (1) P_5 vs sympy Stirling x=5..40, and recurrence symbolically (own transcription from the md)
x=sp.Symbol('x')
P5=x*(x-1)*(x-2)*(x-3)*(x-4)**2*(x-5)**2*(3*x**2-23*x+38)/11520
bad=[n for n in range(5,41) if P5.subs(x,n)!=stirling(n,n-5,kind=2)]
print("P5 vs sympy stirling x=5..40 mismatches:",bad)
# also x=0..4 (n-5 negative index): P5 should be 0
print("P5(0..4):",[P5.subs(x,n) for n in range(5)])
Ps={0:sp.Integer(1),1:x*(x-1)/2,2:x*(x-1)*(x-2)*(3*x-5)/24,3:x*(x-1)*(x-2)**2*(x-3)**2/48,
 4:x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760,5:P5}
for r in range(1,6):
    print("recurrence r=%d diff:"%r, sp.expand(Ps[r].subs(x,x+1)-((x+1-r)*Ps[r-1]+Ps[r])))
    print(" P_r vs stirling n=0..40:", [n for n in range(0,41) if Ps[r].subs(x,n)!=(stirling(n,n-r,kind=2) if n-r>=0 else 0)])
# (2) independent numeric: gaussian-rational series mod s^6
class G:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): return G(s.a+o.a,s.b+o.b)
    def __sub__(s,o): return G(s.a-o.a,s.b-o.b)
    def __mul__(s,o): return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def __truediv__(s,o):
        n=o.a*o.a+o.b*o.b; return s*G(o.a/n,-o.b/n)
N=6
def smul(A,B):
    C=[G() for _ in range(N)]
    for i,u in enumerate(A):
        for j,v in enumerate(B):
            if i+j<N: C[i+j]=C[i+j]+u*v
    return C
def sscale(A,c): return [u*c for u in A]
def sadd(A,B): return [u+v for u,v in zip(A,B)]
one=[G(1)]+[G() for _ in range(N-1)]
ts=[G(Fr(1,2),Fr(1,2)),G(1)]+[G() for _ in range(N-2)]   # t = t*+s
us=[G(Fr(1,2),Fr(-1,2)),G(-1)]+[G() for _ in range(N-2)] # 1-t
t2=smul(ts,ts); u2=smul(us,us); W=sadd(t2,u2)
def spow(A,n):
    R=one
    for _ in range(n): R=smul(R,A)
    return R
def bh(j):
    r=Fr(1)
    for k in range(j): r*=(Fr(1,2)-k)
    return r/factorial(j)
def Sser(d,p):
    q=d-p; tot=[G() for _ in range(N)]
    for j in range(1,d+1):
        Cj=bh(j)*factorial(j)*2**d
        for s in range(max(0,j-q),min(j,p)+1):
            c=Cj*int(stirling(p,s,kind=2))*int(stirling(q,j-s,kind=2))
            if c==0: continue
            term=smul(smul(spow(t2,j-s),spow(u2,s)),spow(W,d-j))
            tot=sadd(tot,sscale(term,G(c)))
    return tot
def Rl(ser,l): return (ser[l]*G(factorial(l)))/ser[0]
d_,m_=sp.symbols('d m'); I=sp.I
N5=sp.sympify(pickle.load(open('N5_proved.pkl','rb'))['N5'],locals={'d':d_,'m':m_,'I':I})
def Pi(l,dd):
    r=1
    for j in range(1,l+1): r*=(2*dd-2*j-1)
    return r
tested=0;bad=[]
for dd in (11,13,17,25):
    for pp in range(0,dd+1):
        if dd-pp>pp: continue   # keys p>=q
        ser=Sser(dd,pp)
        if ser[0].a==0 and ser[0].b==0: bad.append((dd,pp,'S(t*)=0')); continue
        r=Rl(ser,5)
        v=sp.expand(N5.subs({d_:dd,m_:pp}))/(3*Pi(5,dd))
        v=sp.expand(v); pa,pb=sp.re(v),sp.im(v)
        tested+=1
        if sp.Rational(pa)!=r.a or sp.Rational(pb)!=r.b: bad.append((dd,pp))
print("R5 own-series vs N5/(3Pi5): tested",tested,"bad",bad)
# also keys p<q (m < d/2) for d=11,13 to probe 'all (d,m)'
b2=[];t2n=0
for dd in (11,13):
    for pp in range(1,dd//2+1):
        ser=Sser(dd,pp); r=Rl(ser,5)
        v=sp.expand(N5.subs({d_:dd,m_:pp})/(3*Pi(5,dd)))
        t2n+=1
        if sp.Rational(sp.re(v))!=r.a or sp.Rational(sp.im(v))!=r.b: b2.append((dd,pp))
print("p<q keys d=11,13: tested",t2n,"mismatch",b2)
# R1,R2 vs md closed forms
for dd,pp in ((11,6),(17,9)):
    ser=Sser(dd,pp)
    n1=4*pp+(-2+6*I)*dd-4*I*dd**2; v1=sp.expand(n1/(2*dd-3))
    r1=Rl(ser,1)
    print("R1 check",dd,pp, sp.Rational(sp.re(v1))==r1.a and sp.Rational(sp.im(v1))==r1.b)
