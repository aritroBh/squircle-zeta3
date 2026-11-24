import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=120
from fractions import Fraction as F
import engine as E, mzv_exact as M
from mzv_vals import dz
z35=dz(3,5,N=800,J=60); z53=dz(5,3,N=800,J=60)
print('stuffle z35+z53+z8-z3z5:', nstr(z35+z53+zeta(8)-zeta(3)*zeta(5),3))
def monval(mon):
    v=mpf(1)
    for L in mon:
        if len(L)==1: v*=zeta(int(L[0]))
        elif L==(3,5): v*=z35
        else: raise ValueError(L)
    return v
def cnum(a,W):
    r=E.reduced_coeffs(a,8)[W]
    T=E.table(W)
    return sum(mpf(int(c.numerator))/int(c.denominator)*monval(m) for c,m in zip(r,T.monomials))
PILOT='<project>/research/02_ideas/lenses/experimental_pilots'
files={'chain4_coeffs.txt':(1,1,1,1),
       'chain4_generic_1_1_1_2_d105.txt':(1,1,1,2),
       'chain4_generic_1_1_1_2_d75.txt':(1,1,1,2),
       'chain4_generic_1_2_2_1_d75.txt':(1,2,2,1),
       'chain4_generic_2_1_1_2_d75.txt':(2,1,1,2),
       'chain4_generic_3_1_2_5.txt':(3,1,2,5)}
for f,a in files.items():
    p=os.path.join(PILOT,f) if not os.path.exists(f) else f
    if not os.path.exists(p): p=os.path.join(PILOT,f)
    data={}
    for line in open(p):
        parts=line.split()
        if len(parts)>=2 and parts[0].isdigit(): data[int(parts[0])]=mpf(parts[1])
    worst=[]
    for k in range(2,9):
        if k in data:
            d=cnum(a,k)-data[k]
            worst.append((k,nstr(abs(d),2)))
    print(f, a, worst)
print('t=%.1fs'%(time.time()-t0))
