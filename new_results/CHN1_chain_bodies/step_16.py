import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=80
import mzv_exact as M, mzv_num as MN
from mzv_vals import dz
from fractions import Fraction as F
z35=dz(3,5); z53=dz(5,3)
assert abs(z35-mpf('0.20466113696507743533'))<mpf('1e-19')
assert abs(z53-mpf('0.03770767298484754401'))<mpf('1e-19')
print('anchors OK:', nstr(z35,22), nstr(z53,22))
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
print('--- depth-2 W8 words vs dz (target 60+ digits) ---')
for w in [(2,6),(3,5),(4,4),(5,3),(6,2),(7,1)]:
    val=dz(*w); pr=pred8(w)
    print(w, '|diff| =', nstr(abs(val-pr),3))
print('t=%.1fs'%(time.time()-t0))
