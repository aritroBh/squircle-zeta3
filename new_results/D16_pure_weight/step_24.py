# step_23.py -- fixed route-3 (D11c operator formula) + comparison vs D11's exact a_n in d.
import json, sympy as sp
from sympy import Rational, factorial, binomial
from d16_engine import *
base='research/'
KMAX=12
mA={n:{(n,):Rational((-1)**(n+1))/n} for n in range(2,KMAX+1)}
Kexp=s_exp(mA,KMAX)
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
        Bser=bs_add(Bser,{(n-r,r):{(n,):Rational((-1)**n)*2*binomial(n,r)/n}})
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
coshser={0:{():Rational(1)}}
for k in range(1,KMAX//2+1): coshser[2*k]={():Rational(1,factorial(2*k))}
logcosh=s_log(coshser,KMAX)
sqrtcosh=s_exp(s_scale({0:{():Rational(0)},**logcosh},Rational(1,2)),KMAX)
c={k:sqrtcosh.get(2*k,{}).get((),Rational(0)) for k in range(KMAX//2+1)}
print("c_k:",[str(c[k]) for k in range(7)])  # expect 1, 1/4, -1/96, 19/5760, -559/645120...
S={}
for k in range(KMAX//2+1):
    for i in range(0,KMAX+1):
        v=Hser.get((i,2*k),{})
        if v: S=s_add(S,{2*k+i:s_scale({0:v},c[k]*factorial(2*k)/Rational(2)**i)[0]})
S=s_mul(Kexp,S,KMAX)
Pd11c=s_scale(S,4)
Pd11c_tau={k:{m:cc*((-1)**k) for m,cc in d.items()} for k,d in Pd11c.items()}
def from_json(J): return {int(k):{tuple(int(x) for x in m.split()):sp.Rational(cc) for m,cc in d.items()} for k,d in J.items()}
Pmine=from_json(json.load(open('tau_series_2d.json'))['P_over_sqrt2'])
bad=[k for k in range(13) if Pd11c_tau.get(k,{})!=Pmine.get(k,{})]
print("[CHECK 2] D11c operator formula == my B_2 formula (13 coeffs exact):", not bad, bad)
# ---- D11 a_n in d from step1_exact.txt ----
txt=open(base+'08_theory/D11_rhombus_resurgence/step1_exact.txt').read()
import re
seg=txt[txt.find('a_0'):]
print(seg[:1200])
