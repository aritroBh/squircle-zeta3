# step_27.py -- independent numeric engine (mpmath, moment formula as series) to order 64; cross-check vs H1 hi_coeffs.
import mpmath as mp
import json
mp.mp.dps=200
N=64
def ser_exp(A):  # A list, A[0]=0
    R=[mp.mpf(0)]*(N+1); R[0]=mp.mpf(1)
    for n in range(1,N+1):
        R[n]=sum(k*A[k]*R[n-k] for k in range(1,n+1))/n
    return R
def ser_mul(A,B):
    R=[mp.mpf(0)]*(N+1)
    for i in range(N+1):
        if A[i]==0: continue
        for j in range(N+1-i):
            if B[j]!=0: R[i+j]+=A[i]*B[j]
    return R
def logGamma_series(a):
    # series of log Gamma(1 + a*tau) in tau: -gamma*a*tau + sum_{k>=2} (-1)^k zeta(k) a^k tau^k/k
    A=[mp.mpf(0)]*(N+1)
    A[1]=-mp.euler*a
    for k in range(2,N+1):
        A[k]=((-1)**k)*mp.zeta(k)*a**k/k
    return A
def Gamma_series(a):
    return ser_exp(logGamma_series(a))
# B_2(tau) = sum_j C(1/2,j) EX[j], EX[j] = sum_{s<=j} C(j,s)(-1)^{j-s} 2^{-s} sum_{r=0}^s C(s,r) M(r,s)
# M(r,s) = Gamma(1+a+tau(2r-s)) Gamma(1+a-tau(2r-s)) / Gamma(1+a)^2, a=-tau/2
# Note 1+a+tau(2r-s) = 1 + tau(2r-s-1/2): linear in tau => Gamma series with coefficient (2r-s-1/2)
from math import comb, factorial
def halfbinom(j):
    v=mp.mpf(1)
    for i in range(j): v*=mp.mpf(1)/2-i
    return v/factorial(j)
B=[mp.mpf(0)]*(N+1); B[0]=mp.mpf(1)
for j in range(1,N+1):
    cj=halfbinom(j)
    EX=[mp.mpf(0)]*(N+1)
    for s in range(0,j+1):
        RS=[mp.mpf(0)]*(N+1)
        for r in range(0,s+1):
            M=ser_mul(Gamma_series(mp.mpf(2*r-s)-mp.mpf(1)/2), Gamma_series(-(mp.mpf(2*r-s)-mp.mpf(1)/2)))
            for k in range(N+1): RS[k]+=comb(s,r)*M[k]
        Ginv_a=Gamma_series(-mp.mpf(1)/2)   # Gamma(1+a) = Gamma(1-tau/2)
        Ginv_a2=ser_mul(Ginv_a,Ginv_a)
        # divide RS by Ginv_a2: series division
        Q=[mp.mpf(0)]*(N+1)
        for n_ in range(N+1):
            Q[n_]=(RS[n_]-sum(Ginv_a2[k]*Q[n_-k] for k in range(1,n_+1)))/Ginv_a2[0]
        fac=(-1)**(j-s)*comb(j,s)*mp.mpf(2)**(-s)
        for k in range(N+1): EX[k]+=fac*Q[k]
    for k in range(N+1): B[k]+=cj*EX[k]
    if j%16==0: print("j=",j,"done")
# S_2/S_2(1) = Gamma(1-tau/2)^2 / Gamma(1-tau) * B
G1=Gamma_series(-mp.mpf(1)/2); G2=Gamma_series(mp.mpf(-1))
num=ser_mul(ser_mul(G1,G1),B)
S=[mp.mpf(0)]*(N+1)
for n_ in range(N+1):
    S[n_]=(num[n_]-sum(G2[k]*S[n_-k] for k in range(1,n_+1)))/G2[0]
json.dump([mp.nstr(x,180) for x in S], open('hi_tau_numeric_n2.json','w'))
print("S computed. S[2..6]:",[mp.nstr(x,40) for x in S[2:7]])
# compare with H1
hj=json.load(open('research/10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json'))
worst=mp.mpf(0); wk=0
for k in range(2,65):
    hv=mp.mpf(hj['eta'][k])
    rel=abs((S[k]-hv)/hv)
    if rel>worst: worst=rel; wk=k
print("worst relative deviation vs H1 hi_coeffs (k=2..64):", mp.nstr(worst,5), "at k=",wk)
