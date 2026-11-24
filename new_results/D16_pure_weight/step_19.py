# step_18.py -- exact engine v3: pure rational dict-arithmetic over QQ[Z_2..Z_K] truncated by weight.
# A series is dict: tau-power -> dict: z-monomial (tuple of z-indices, sorted) -> Rational.
# Weight of monomial (i1,...,ir) = i1+...+ir ; only terms with tau-power == total weight needed,
# but we keep general (power, monomial) since B_n has pure weight anyway (checked).
import sympy as sp
from sympy import Rational, factorial
from sympy.utilities.iterables import partitions as spart
import time, json

KMAX = 12
ONE = {0: {(): Rational(1)}}

def s_add(A,B):
    R={}
    for src in (A,B):
        for p,md in src.items():
            d=R.setdefault(p,{})
            for m,c in md.items(): d[m]=d.get(m,0)+c
    return {p:{m:c for m,c in d.items() if c} for p,d in R.items() if any(d.values())}
def s_scale(A,q):
    return {p:{m:c*q for m,c in d.items()} for p,d in A.items()}
def s_mul(A,B,K):
    R={}
    for p1,d1 in A.items():
        for p2,d2 in B.items():
            if p1+p2>K: continue
            d=R.setdefault(p1+p2,{})
            for m1,c1 in d1.items():
                for m2,c2 in d2.items():
                    m=tuple(sorted(m1+m2))
                    d[m]=d.get(m,0)+c1*c2
    return {p:{m:c for m,c in d.items() if c} for p,d in R.items() if any(d.values())}
def s_exp(A,K):  # A has no constant term
    R={0:{():Rational(1)}}
    for n in range(1,K+1):
        d={}
        for k in range(1,n+1):
            if k in A and n-k in R:
                for m,c in A[k].items():
                    for m2,c2 in R[n-k].items():
                        mm=tuple(sorted(m+m2)); d[mm]=d.get(mm,0)+Rational(k,n)*c*c2
        if d: R[n]={m:c for m,c in d.items() if c}
    return R
def s_log(A,K):  # A[0]=1
    L={}
    for n in range(1,K+1):
        d={}
        if n in A:
            for m,c in A[n].items(): d[m]=d.get(m,0)+n*c
        for k in range(1,n):
            if k in L and n-k in A:
                for m,c in L[k].items():
                    for m2,c2 in A[n-k].items():
                        mm=tuple(sorted(m+m2)); d[mm]=d.get(mm,0)-k*c*c2
        dd={m:sp.nsimplify(0)+Rational(c)/n for m,c in d.items() if c}
        L[n]={m:Rational(c)/n for m,c in d.items() if c}
    return L

def zterm(k, q):
    return {k:{(k,):Rational(q)}}

def build_CS(n, K):
    """CS(s): dict s -> series. Exact rational arithmetic; n is a python int."""
    CS={0:{0:{():Rational(1)}}}
    for s in range(1,K+1):
        acc={}
        for p in spart(s):
            parts=[]; mult=[]
            for part,m in sorted(p.items()):
                parts += [part]*m; mult.append(m)
            r=len(parts)
            coef = Rational(factorial(s))
            for t in range(r): coef *= (n-t)
            for m in mult: coef /= factorial(m)
            for l in parts: coef /= factorial(l)
            cvals = [Rational(2*l*n - 2*s - (n-1), n) for l in parts]
            c0 = Rational(-2*s-(n-1), n)
            cm = Rational(-(n-1), n)
            lm={}
            for k in range(2,K+1):
                br = sum(c**k for c in cvals) + (n-r)*c0**k - n*cm**k
                if br: lm[k]={(k,):Rational((-1)**k)*br/k}
            M=s_exp(lm,K)
            acc=s_add(acc, s_scale(M,coef))
        CS[s]=acc
    return CS

def halfbinom(j):
    v=Rational(1)
    for i in range(j): v*=Rational(1,2)-i
    return v/factorial(j)

def build_B(n,K,CS):
    B={0:{():Rational(1)}}
    for j in range(1,K+1):
        cj=halfbinom(j)
        EX={}
        for s in range(0,j+1):
            fac = Rational((-1)**(j-s))*sp.binomial(j,s)*Rational(n)**(j-s)/Rational(n)**j
            EX=s_add(EX, s_scale(CS[s],fac))
        B=s_add(B, s_scale(EX,cj))
    return B

t0=time.time()
CS2=build_CS(2,KMAX); B2=build_B(2,KMAX,CS2)
print("exact B_2, time",round(time.time()-t0,1),"s")
def show(B,name):
    for k in sorted(B):
        if k==0: continue
        terms=" + ".join(f"({c})*{'*'.join('Z%d'%i for i in m) if m else '1'}" for m,c in sorted(B[k].items()))
        print(f"[tau^{k}] {name} = {terms}")
show(B2,'B_2')
# quick check vs known: [tau^2] should be Z2/2, [tau^3] Z3/2, [tau^4] -Z2^2/8 + Z4/4
assert B2[2]=={(2,):Rational(1,2)}, B2[2]
assert B2[3]=={(3,):Rational(1,2)}, B2[3]
assert B2[4]=={(2,2):Rational(-1,8),(4,):Rational(1,4)}, B2[4]
assert B2[5]=={(2,3):Rational(-1,4)}, B2[5]
print("low-order assertions OK")
