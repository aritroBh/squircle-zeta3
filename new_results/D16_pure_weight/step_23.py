# step_22.py -- route 3: D11c operator formula P/sqrt2 = 4*K*T (eta_D = -tau), and check vs D11's a_n in d.
import json, sympy as sp
from sympy import Rational, factorial, binomial
from d16_engine import *
KMAX=12
# K(eta_D) = exp(-sum (-1)^n Z_n eta_D^n / n), eta_D = 1/p - 1 = -tau
mA={n:{(n,):Rational((-1)**(n+1))/n} for n in range(2,KMAX+1)}     # in eta_D
Kexp=s_exp(mA,KMAX)
# H(mu,s) = exp( sum_n (-1)^n Z_n [(mu+s)^n + (mu-s)^n]/n ) ; need d_s^{2k} H at s=0, mu=eta_D/2
# work in biseries: powers (i,j) of (mu,s)
def bs_add(A,B):
    R={}
    for src in (A,B):
        for key,d in src.items():
            dd=R.setdefault(key,{})
            for m,c in d.items(): dd[m]=dd.get(m,0)+c
    return R
Bser={}
for n in range(2,KMAX+1):
    for r in range(0,n+1,2):
        cf=Rational((-1)**n)*2*binomial(n,r)/n
        Bser=bs_add(Bser,{(n-r,r):{(n,):cf}})
def bs_exp(A,K):
    R={(0,0):{():Rational(1)}}
    for tot in range(1,K+1):
        for i in range(0,tot+1):
            j=tot-i; d={}
            if i>0:
                for (i1,j1),aa in A.items():
                    if i1>=1 and (i-i1,j-j1) in R:
                        for m,c in aa.items():
                            for m2,c2 in R[(i-i1,j-j1)].items():
                                mm=tuple(sorted(m+m2)); d[mm]=d.get(mm,0)+Rational(i1,i)*c*c2
            else:
                for (i1,j1),aa in A.items():
                    if i1==0 and j1>=1 and (0,j-j1) in R:
                        for m,c in aa.items():
                            for m2,c2 in R[(0,j-j1)].items():
                                mm=tuple(sorted(m+m2)); d[mm]=d.get(mm,0)+Rational(j1,j)*c*c2
            if d: R[(i,j)]={m:c for m,c in d.items() if c}
    return R
Hser=bs_exp(Bser,KMAX)
# c_k: sqrt(cosh z)
coshs={2*k:{():Rational(1,factorial(2*k))} for k in range(1,KMAX//2+1)}
# sqrt via log/exp
cl=s_log({0:{():Rational(1)},**coshs},KMAX)
sqc=s_exp({k:s_scale({0:{():Rational(0)},k:v},Rational(1,2))[k] for k,v in cl.items()},KMAX)
c={k:sqc.get(2*k,{}).get((),Rational(0)) for k in range(KMAX//2+1)}
# T(eta_D) = sum_k c_k eta_D^{2k} (d_s)^{2k} H (mu=eta_D/2) ; S = K*T ; P/sqrt2 = 4 S
S={0:{():Rational(1)}}
for k in range(KMAX//2+1):
    for i in range(0,KMAX+1):
        v=Hser.get((i,2*k),{})
        if v:
            fac=c[k]*factorial(2*k)/Rational(2)**i
            S=s_add(S,{2*k+i:s_scale({0:v},fac)[0]})
S=s_mul(Kexp,S,KMAX)
Pd11c=s_scale(S,4)   # P/sqrt2 as series in eta_D
# convert to tau: eta_D = -tau
Pd11c_tau={k:{m:c*((-1)**k) for m,c in d.items()} for k,d in Pd11c.items()}
mine=json.load(open('tau_series_2d.json'))['P_over_sqrt2']
def from_json(J):
    return {int(k):{tuple(int(x) for x in m.split()):sp.Rational(c) for m,c in d.items()} for k,d in J.items()}
Pmine=from_json(mine)
ok=all((Pd11c_tau.get(k,{})==Pmine.get(k,{})) for k in range(13))
print("[CHECK 2] D11c operator formula (route 3) == my B_2 formula (route 1):",ok)
for k in range(13):
    if Pd11c_tau.get(k,{})!=Pmine.get(k,{}): print("k=",k,Pd11c_tau.get(k),Pmine.get(k))
# ---- check vs D11's a_n (d-series): substitute tau = d/(1+d) ----
# tau^j = d^j (1+d)^{-j} = sum_{m>=j} (-1)^{m-j} C(m-1,j-1) d^m
import itertools
dser={0:{():Rational(4)}}
for j in range(2,KMAX+1):
    for m in range(j,KMAX+1):
        cf=Rational((-1)**(m-j))*sp.binomial(m-1,j-1)
        dser=s_add(dser,{m:s_scale({0:Pmine[j]},cf)[0]})
# D11 a_n (reduced, zetas): from D11 result.md / step1_exact.txt
import os
f=open(base+'08_theory/D11_rhombus_resurgence/step1_exact.txt').read()
i0=f.find('a_2'); print("\nD11 step1_exact.txt snippet:\n",f[f.find('a_0'):f.find('a_0')+900])
