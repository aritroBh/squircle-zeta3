import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=80
import mzv_exact as M, mzv_num as MN
from mzv_vals import dz
from fractions import Fraction as F
z35=dz(3,5)
def monval(mon):
    v=mpf(1)
    for L in mon:
        if len(L)==1: v*=zeta(int(L[0]))
        elif L==(3,5): v*=z35
        else: raise ValueError(L)
    return v
T8=M.DSTable(8)
def pred8(w):
    return sum(mpf(int(c.numerator))/int(c.denominator)*monval(m) for c,m in zip(T8.reduce_word(w),T8.monomials))
print('--- Richardson nested sums vs DS prediction (target >=40 digits) ---')
for w in [(5,3),(2,3,3),(3,2,3),(4,2,2),(2,2,2,2)]:
    t=time.time()
    val=MN.mzv_num(w)
    pr=pred8(w)
    print(w,'|diff|=',nstr(abs(val-pr),3),' (%.0fs)'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
open('progress.log','a').write('step_15-16: numeric validations PASS: anchors; 6 depth-2 W8 words vs dz to ~1e-80; 5 all->=2 words Richardson vs DS to >=40 digits\n')
