# step_28.py -- fast independent numeric engine: cumulant route. B_2 = sum_k c_k tau^{2k} E[l^{2k}](tau).
import mpmath as mp, json
from math import factorial, comb
mp.mp.dps=150
N=64
def ser_exp(A):
    R=[mp.mpf(0)]*(N+1); R[0]=mp.mpf(1)
    for n in range(1,N+1): R[n]=sum(k*A[k]*R[n-k] for k in range(1,n+1))/n
    return R
def ser_mul(A,B):
    R=[mp.mpf(0)]*(N+1)
    for i in range(N+1):
        ai=A[i]
        if ai==0: continue
        for j in range(N+1-i):
            if B[j]: R[i+j]+=ai*B[j]
    return R
# cumulants of l = log(G1/G2), G~Gamma(1-tau/2): kap_{2m}(tau) = 2 sum_{t>=0} tau^t (2m-1+t)! zeta(2m+t)/(2^t t!)
def kap(m):
    A=[mp.mpf(0)]*(N+1)
    for t in range(0,N+1):
        A[t]=2*mp.mpf(factorial(2*m-1+t))*mp.zeta(2*m+t)/(mp.mpf(2)**t*factorial(t))
    return A
# moments E[l^{2k}](tau) via exp of the cgf in s
El={0:[mp.mpf(1)]+[mp.mpf(0)]*N}
ksers={m:kap(m) for m in range(1,N//2+1)}
# exp: E[l^{2k}]/(2k)! = coeff of s^{2k} in exp(sum kap_{2m} s^{2m}/(2m)!)
# compute iteratively: M[2k] = sum_{m<=k} kap_{2m} * binom(2k-1,2m-1) * M[2k-2m]
El[0]=[mp.mpf(1)]+[mp.mpf(0)]*N
for k in range(1,N//2+1):
    acc=[mp.mpf(0)]*(N+1)
    for m in range(1,k+1):
        cm=ksers[m]
        w=comb(2*k-1,2*m-1)
        prev=El[k-m]
        for t in range(N+1): acc[t]+=w*cm[t]*prev[0]*0 or acc[t]  # placeholder
    El[k]=acc
# do it properly with series convolution:
El={0:[mp.mpf(1)]+[mp.mpf(0)]*N}
for k in range(1,N//2+1):
    acc=[mp.mpf(0)]*(N+1)
    for m in range(1,k+1):
        cm=ksers[m]; w=comb(2*k-1,2*m-1); prev=El[k-m]
        for t in range(N+1):
            cmt=cm[t]
            if cmt==0: continue
            for u in range(N+1-t):
                if prev[u]: acc[t+u]+=w*cmt*prev[u]
    El[k]=acc
# sqrt(cosh) coefficients c_k (rational)
c=[]
import sympy as sp
z=sp.symbols('z'); sersp=sp.series(sp.sqrt(sp.cosh(z)),z,0,N+2).removeO()
c=[sp.Rational(sersp.coeff(z,2*k)) for k in range(N//2+1)]
# B_2(tau) = sum_k c_k tau^{2k} El[k](tau)
B=[mp.mpf(0)]*(N+1); B[0]=mp.mpf(1)
for k in range(1,N//2+1):
    ck=mp.mpf(str(c[k]))
    for t in range(N+1-2*k): B[t+2*k]+=ck*El[k][t]
# S = Gamma(1-tau/2)^2/Gamma(1-tau) * B
def logGamma_series(a):
    A=[mp.mpf(0)]*(N+1); A[1]=-mp.euler*a
    for k in range(2,N+1): A[k]=((-1)**k)*mp.zeta(k)*a**k/k
    return A
G1=ser_exp(logGamma_series(-mp.mpf(1)/2)); G2=ser_exp(logGamma_series(mp.mpf(-1)))
num=ser_mul(ser_mul(G1,G1),B)
S=[mp.mpf(0)]*(N+1)
for n_ in range(N+1): S[n_]=(num[n_]-sum(G2[k]*S[n_-k] for k in range(1,n_+1)))/G2[0]
json.dump([mp.nstr(x,130) for x in S], open('hi_tau_numeric_n2.json','w'))
hj=json.load(open('research/10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json'))
worst=mp.mpf(0); wk=0
for k in range(2,65):
    hv=mp.mpf(hj['eta'][k]); rel=abs((S[k]-hv)/hv)
    if rel>worst: worst, wk = rel, k
print("worst rel deviation vs H1 hi_coeffs_n2 (k=2..64):", mp.nstr(worst,5), "at k=",wk)
print("S[64]:", mp.nstr(S[64],30))
print("H1   :", hj['eta'][64][:40])
