
import sympy as sp
from sympy import Rational, factorial
from sympy.utilities.iterables import partitions as spart

def s_add(A,B):
    R={}
    for src in (A,B):
        for p,md in src.items():
            d=R.setdefault(p,{})
            for m,c in md.items(): d[m]=d.get(m,0)+c
    return {p:{m:c for m,c in d.items() if c} for p,d in R.items() if any(d.values())}
def s_scale(A,q):
    return {p:{m:Rational(c)*q for m,c in d.items()} for p,d in A.items()}
def s_mul(A,B,K):
    R={}
    for p1,d1 in A.items():
        for p2,d2 in B.items():
            if p1+p2>K: continue
            d=R.setdefault(p1+p2,{})
            for m1,c1 in d1.items():
                for m2,c2 in d2.items():
                    m=tuple(sorted(m1+m2)); d[m]=d.get(m,0)+c1*c2
    return {p:{m:c for m,c in d.items() if c} for p,d in R.items() if any(d.values())}
def s_exp(A,K):
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
def s_log(A,K):
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
        L[n]={m:Rational(c)/n for m,c in d.items() if c}
    return L
def halfbinom(j):
    v=Rational(1)
    for i in range(j): v*=Rational(1,2)-i
    return v/factorial(j)
def build_CS(n,K):
    """composition sums CS(s), s=0..K ; n python int >= 2. Exact."""
    CS={0:{0:{():Rational(1)}}}
    for s in range(1,K+1):
        acc={}
        for p in spart(s):
            parts=[]; mult=[]
            for part,m in sorted(p.items()):
                parts += [part]*m; mult.append(m)
            r=len(parts)
            if r>n: continue
            coef = Rational(factorial(s))
            for t in range(r): coef *= (n-t)
            for m in mult: coef /= factorial(m)
            for l in parts: coef /= factorial(l)
            cvals=[Rational(2*l*n-2*s-(n-1),n) for l in parts]
            c0=Rational(-2*s-(n-1),n); cm=Rational(-(n-1),n)
            lm={}
            for k in range(2,K+1):
                br=sum(c**k for c in cvals)+(n-r)*c0**k-n*cm**k
                if br: lm[k]={(k,):Rational((-1)**k)*br/k}
            M=s_exp(lm,K)
            acc=s_add(acc,s_scale(M,coef))
        CS[s]=acc
    return CS
def build_B(n,K,CS=None):
    if CS is None: CS=build_CS(n,K)
    B={0:{():Rational(1)}}
    for j in range(1,K+1):
        cj=halfbinom(j)
        EX={}
        for s in range(0,j+1):
            fac=Rational((-1)**(j-s))*sp.binomial(j,s)*Rational(n)**(j-s)/Rational(n)**j
            EX=s_add(EX,s_scale(CS[s],fac))
        B=s_add(B,s_scale(EX,cj))
    return {p:d for p,d in B.items() if d}
def Lambda(n,K):
    """Lambda_n(tau) = -sum_{k>=2} Z_k (n-1)^k (1-n^{1-k}) tau^k/k"""
    return {k:{(k,):-Rational((n-1)**k)*(1-Rational(1,n**(k-1)))/k} for k in range(2,K+1)}
def logRat(n,K):
    """log Rat_n(tau) = (n-1) log(1-tau) - sum_{i=1}^{n-1} log(1-(n-1)tau/i), rational coeffs"""
    L={}
    for k in range(1,K+1):
        L[k]={(): -(n-1)/Rational(k) + sum(Rational((n-1)**k)/(Rational(i)**k*k) for i in range(1,n))}
    return L
def full_S(n,K):
    """returns S_n/S_n(1), log(S_n/(S_n(1) Rat_n)) = Lambda + log B, and parts"""
    CS=build_CS(n,K); B=build_B(n,K,CS)
    eLam=s_exp(Lambda(n,K),K)
    S=s_mul(eLam,B,K)  # = S_n/S_n(1)/Rat_n ... wait: S_n/S_n(1) = Rat_n * e^Lambda * B
    return B,eLam,S
def show(S,name,powers=None):
    for k in sorted(S):
        if k==0: continue
        terms=" + ".join(f"({c})*"+("*".join('Z%d'%i for i in m) if m else '1') for m,c in sorted(S[k].items()))
        print(f"[tau^{k}] {name} = {terms}")
def ser_to_json(S): return {str(k):{(' '.join(map(str,m))):str(c) for m,c in d.items()} for k,d in S.items()}
