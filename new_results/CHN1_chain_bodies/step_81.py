import os, time, json
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=90
from fractions import Fraction as F
import mzv_exact as M, engine as E
from mzv_vals import dz
def table(W):
    if W not in E._TABLES: E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
z39=dz(3,9,N=800,J=60); z57=dz(5,7,N=800,J=60); z35=dz(3,5,N=800,J=60); z37=dz(3,7,N=800,J=60)
r=E.reduced_coeffs((1,1,1,1),12)[12]
def monval(mon):
    v=mpf(1)
    for L in mon:
        if len(L)==1: v*=zeta(int(L[0]))
        elif L==(3,9): v*=z39
        elif L==(5,7): v*=z57
        elif L==(3,5): v*=z35
        elif L==(3,7): v*=z37
        else: raise ValueError(L)
    return v
T=E.table(12)
pred=sum(mpf(int(c.numerator))/int(c.denominator)*monval(m) for c,m in zip(r,T.monomials) if c!=0)
open('diag_c12_prediction.txt','w').write(nstr(pred,80)+'\n')
print('predicted diag c12 =', nstr(pred,40))
print('job done?', os.path.exists('diag_c12_result.json'))
if os.path.exists('diag_c12_result.json'):
    res=json.load(open('diag_c12_result.json'))
    print('cauchy c12 =', res['12'][0][:55], ' imag', res['12'][1])
    d=mpf(res['12'][0])-pred
    print('|diff| =', nstr(abs(d),3))
print('t=%.1fs'%(time.time()-t0))
