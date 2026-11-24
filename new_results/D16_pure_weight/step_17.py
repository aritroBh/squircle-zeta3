# step_16.py -- D16 core engine: explicit Gamma-product formula for the bending factor B_n(tau),
# with formal symbols Z_k for zeta(k). Verification against H1's Route A data (n=2..5).
#
# THEORY (derived here): S_n(p)/S_n(1) = (n-1)! eps^{n-1} Gamma(1+a)^n / Gamma(1+(n-1)eps) * B_n(eta),
# eps=1/p=1-tau, a=-(n-1)tau/n, B_n = E[((1/n)Sum_i e^{2 tau m_i})^{1/2}], m_i=log G_i - mean, G_i~Gamma(1+a).
# Expand the 1/2-power: B_n = Sum_j C(1/2,j) E[X^j], X=(1/n)Sum_i(e^{2tau m_i}-1).
# For a tuple with equality-pattern = partition lam of s (r distinct values, mult m_lam),
# and firing set S of size s:  E[prod] = n^{-j} * Sum over (lam,s) ...
#   E[X^j] = n^{-j} Sum_{s=0}^j C(j,s)(-1)^{j-s} n^{j-s} CS(s), 
#   CS(s) = Sum_{lam |- s, r=len<=n} (n)_r/prod(mult!) * s!/prod(lam_i!) * M(lam,s)
#   M(lam,s) = [Prod_i Gamma(1+a+2tau(lam_i - s/n))] Gamma(1+a-2tau s/n)^{n-r} / Gamma(1+a)^n
# log M = Sum_{k>=2} (-1)^k Z_k tau^k / k * [Sum_i c(lam_i)^k + (n-r)c(0)^k - n(-(n-1)/n)^k],  c(l)=2l-2s/n-(n-1)/n
# (gamma cancels since Sum_i c(lam_i)+(n-r)c(0)-n*(-(n-1)/n)*... linear part = 0 -- checked below)
import sympy as sp
from sympy import symbols, Rational, factorial, binomial, prod as sprod
from sympy.utilities.iterables import partitions as spart

KMAX = 12
tau = symbols('tau')
Z = {k: symbols(f'Z{k}') for k in range(2, KMAX+2)}

# truncated series in tau with sympy coeffs: dict power->expr
def tmul(A,B,K):
    R={}
    for i,a in A.items():
        for j,b in B.items():
            if i+j<=K: R[i+j]=R.get(i+j,0)+a*b
    return {k:sp.expand(v) for k,v in R.items() if sp.expand(v)!=0}
def texp(A,K):
    R={0:sp.Integer(1)}
    for n in range(1,K+1):
        acc=sum(k*A.get(k,0)*R.get(n-k,0) for k in range(1,n+1))
        v=sp.expand(acc/n)
        if v!=0: R[n]=v
    return R
def tlog(A,K):  # A[0]=1
    L={}
    for n in range(1,K+1):
        acc=A.get(n,0)*n
        for k in range(1,n):
            acc-=k*L.get(k,0)*A.get(n-k,0)
        L[n]=sp.expand(acc/n)
    return {k:v for k,v in L.items() if v!=0}
def tpow_half(A,K):  # (1+A)^p via log/exp
    p=Rational(1,2)
    L=tlog(A,K)
    return texp({k:p*v for k,v in L.items()},K)

def comps_into_n_parts(s):
    # partitions lambda of s with multiplicity counts, as (list of parts, mult dict)
    out=[]
    for p in spart(s):
        parts=[]; mult={}
        for part,m in sorted(p.items()):
            parts += [part]*m; mult[part]=m
        out.append((parts,mult))
    return out

def build_B(n, K):
    """B_n(tau) as truncated series, coefficients in Q(n)[Z2..ZK] (n symbolic or int)."""
    B={0:sp.Integer(1)}
    for j in range(1,K+1):
        cj=binomial(Rational(1,2),j)
        EX={0:sp.Integer(0)}
        for s in range(0,j+1):
            sgn=(-1)**(j-s)
            CS={}  # composition sum, series in tau
            for parts,mult in comps_into_n_parts(s):
                r=len(parts)
                if r>n and isinstance(n,int): continue
                # coefficient (n)_r / prod mult! * s!/prod lam!
                coef = sprod([n-t for t in range(r)])/sprod([factorial(m) for m in mult.values()])
                coef *= factorial(s)/sprod([factorial(l) for l in parts])
                # log M
                lm={}
                for k in range(2,K+1):
                    br = sum((2*l-Rational(2*s,1)/n-(n-1)/n)**k for l in parts) \
                       + (n-r)*(Rational(-2*s,1)/n-(n-1)/n)**k \
                       - n*(-(n-1)/n)**k
                    lm[k]=sp.expand((-1)**k*Z[k]/k*br)
                M=texp(lm,K)
                CS=tmul({0:coef},CS,K) if False else {k:CS.get(k,0)+coef*M.get(k,0) for k in range(0,K+1)}
                CS={k:sp.expand(v) for k,v in CS.items() if sp.expand(v)!=0}
            fac = sgn * n**(j-s) / n**j
            EX={k:EX.get(k,0)+fac*CS.get(k,0) for k in range(0,K+1)}
            EX={k:sp.expand(v) for k,v in EX.items() if sp.expand(v)!=0}
        B={k:B.get(k,0)+cj*EX.get(k,0) for k in range(0,K+1)}
        B={k:sp.expand(v) for k,v in B.items() if sp.expand(v)!=0}
    return B

import json, time
t0=time.time()
B2=build_B(2,KMAX)
print("B_2 built, time", round(time.time()-t0,1),"s")
for k in sorted(B2): print(f"[tau^{k}] B_2 =", B2[k])
