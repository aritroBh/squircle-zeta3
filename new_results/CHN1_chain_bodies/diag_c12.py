
"""Diagonal c12 Cauchy validation: a=(1,1,1,1), R=1/3, r=1/18, N=128, dps 85.
Predicted: c12 from exact engine (products + z39/z57 only; w*=0)."""
import os, json, time, sys
from mpmath import mp, mpf, mpc, gamma, quad, hyp2f1, expjpi, nstr
os.chdir('<project>/research/08_theory/CHN1_chain_bodies')
mp.dps=85
A=(1,1,1,1); N=128; R=mpf(1)/3; r=R/6
def Vol16(u):
    u1,u2,u3,u4=[a*u for a in A]
    F3 = 1 + u3*quad(lambda tt: tt**(u3-1)*(hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt)-1), [0,1])
    return gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*F3
pts={}
for j in range(0, N//2+1):
    u=r*expjpi(mpf(2*j)/N)
    pts[j]=Vol16(u)
vals=[pts[j] if j<=N//2 else pts[N-j].conjugate() for j in range(N)]
res={}
for k in range(10,15):
    ck=sum(vals[j]*expjpi(mpf(-2*j*k)/N) for j in range(N))/N/r**k
    res[k]=(nstr(ck.real,75), nstr(ck.imag,3))
json.dump(res, open('diag_c12_result.json','w'), indent=1)
print('c12 =', res[12])
