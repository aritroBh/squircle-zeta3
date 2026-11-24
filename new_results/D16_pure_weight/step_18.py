# step_17.py -- exact engine v2: exact half-binomials, cached composition sums CS(s); B_n(tau) with formal Z_k.
import sympy as sp
from sympy import symbols, Rational, factorial, prod as sprod
from sympy.utilities.iterables import partitions as spart
import json, time

KMAX = 12
tau = symbols('tau')
Z = {k: symbols(f'Z{k}') for k in range(2, KMAX+2)}

def halfbinom(j):
    v = Rational(1,1)
    for i in range(j): v *= Rational(1,2)-i
    return v/factorial(j)

def texp(A,K):
    R={0:sp.Integer(1)}
    for n in range(1,K+1):
        v=sp.expand(sum(k*A.get(k,0)*R.get(n-k,0) for k in range(1,n+1))/n)
        if v!=0: R[n]=v
    return R
def tlog(A,K):
    L={}
    for n in range(1,K+1):
        acc=A.get(n,0)*n - sum(k*L.get(k,0)*A.get(n-k,0) for k in range(1,n))
        L[n]=sp.expand(acc/n)
    return {k:v for k,v in L.items() if v!=0}

def build_CS(n, K):
    """CS(s) for s=0..K: composition sums as series in tau with Q(n)[Z] coeffs."""
    CS={}
    for s in range(0,K+1):
        acc={0:sp.Integer(0)}
        if s==0:
            CS[0]={0:sp.Integer(1)}; continue
        for p in spart(s):
            parts=[]; mult=[]
            for part,m in sorted(p.items()):
                parts += [part]*m; mult.append(m)
            r=len(parts)
            coef = sprod([n-t for t in range(r)])/sprod([factorial(m) for m in mult])
            coef *= factorial(s)/sprod([factorial(l) for l in parts])
            lm={}
            for k in range(2,K+1):
                br = sum((2*l-Rational(2*s)/n-(n-1)/n)**k for l in parts) \
                   + (n-r)*(-Rational(2*s)/n-(n-1)/n)**k - n*(-(n-1)/n)**k
                lm[k]=sp.expand((-1)**k*Z[k]*br/k)
            M=texp(lm,K)
            for kk,v in M.items():
                acc[kk]=sp.expand(acc.get(kk,0)+coef*v)
        CS[s]={k:v for k,v in acc.items() if v!=0}
    return CS

def build_B(n,K,CS=None):
    if CS is None: CS=build_CS(n,K)
    B={0:sp.Integer(1)}
    for j in range(1,K+1):
        cj=halfbinom(j)
        EX={}
        for s in range(0,j+1):
            fac = (-1)**(j-s)*sp.binomial(j,s)*n**(j-s)/n**j
            for kk,v in CS[s].items():
                EX[kk]=sp.expand(EX.get(kk,0)+fac*v)
        for kk,v in EX.items():
            if v!=0: B[kk]=sp.expand(B.get(kk,0)+cj*v)
    return {k:v for k,v in B.items() if v!=0}

t0=time.time()
CS2=build_CS(2,KMAX)
B2=build_B(2,KMAX,CS2)
print("B_2 exact, time",round(time.time()-t0,1),"s ; B_2[0] =",B2.get(0))
for k in sorted(B2): print(f"[tau^{k}] B_2 =",B2[k])
json.dump({str(k):str(v) for k,v in B2.items()}, open('B2_exact.json','w'))
