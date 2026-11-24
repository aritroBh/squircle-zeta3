import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=90
from fractions import Fraction as F
import mzv_exact as M, engine as E, mzv_num as MN
from mzv_vals import dz
def table(W):
    if W not in E._TABLES: E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
# numeric basis values W9..W12
z35=dz(3,5,N=800,J=60); z37=dz(3,7,N=800,J=60); z39=dz(3,9,N=800,J=60); z57=dz(5,7,N=800,J=60)
t=time.time(); z335=MN.mzv_num((3,3,5), npts=16, N0=300, step=120); print('z335 eval %.0fs'%(time.time()-t))
t=time.time(); z2118=MN.mzv_num((2,1,1,8)) if False else None
def monval(mon):
    v=mpf(1)
    for L in mon:
        if len(L)==1: v*=zeta(int(L[0]))
        elif L==(3,5): v*=z35
        elif L==(3,7): v*=z37
        elif L==(3,9): v*=z39
        elif L==(5,7): v*=z57
        elif L==(3,3,5): v*=z335
        else: raise ValueError(L)
    return v
def cnum(a,W):
    r=E.reduced_coeffs(a,12)[W]
    T=E.table(W)
    return sum(mpf(int(c.numerator))/int(c.denominator)*monval(m) for c,m in zip(r,T.monomials) if not (m==((2,1,1,8),)))
# NOTE: w* coefficient excluded numerically (no value); on planes/diag w*=0 anyway; for generic pts only report w*=0 checks separately
PILOT='<project>/research/02_ideas/lenses/experimental_pilots'
files={'chain4_coeffs.txt':(1,1,1,1),
       'chain4_generic_1_1_1_2_d105.txt':(1,1,1,2),
       'chain4_generic_1_1_1_2_d75.txt':(1,1,1,2),
       'chain4_generic_1_2_2_1_d75.txt':(1,2,2,1),
       'chain4_generic_2_1_1_2_d75.txt':(2,1,1,2),
       'chain4_generic_3_1_2_5.txt':(3,1,2,5)}
for f,a in files.items():
    p=f if os.path.exists(f) else os.path.join(PILOT,f)
    data={}
    for line in open(p):
        parts=line.split()
        if len(parts)>=2 and parts[0].isdigit(): data[int(parts[0])]=mpf(parts[1])
    diffs=[]
    for k in (9,10,11):
        if k in data:
            # w* coefficient is 0 at W=12 only; for W=9,10,11 all basis elements have values
            diffs.append((k,nstr(abs(cnum(a,k)-data[k]),2)))
    print(f, a, diffs)
print('t=%.1fs'%(time.time()-t0))
