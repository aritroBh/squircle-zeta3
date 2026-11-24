import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr
mp.dps=80
import mzv_exact as M, mzv_num as MN
from mzv_vals import dz
from fractions import Fraction as F
# anchors
z35=dz(3,5); z53=dz(5,3)
print('z(3,5)=',nstr(z35,30)); print('z(5,3)=',nstr(z53,30))
assert nstr(z35,20)=='0.20466113696507743533', nstr(z35,20)
assert nstr(z53,20)=='0.03770767298484754401', nstr(z53,20)
print('anchors OK')
# numeric values of basis monomials
BASVAL={}
def monval(mon):
    v=mpf(1)
    for L in mon:
        if len(L)==1: v*=zeta(int(L[0]))
        elif L==(3,5): v*=z35
        elif L==(3,7): v*=dz(3,7)
        else: raise ValueError(L)
    return v
T8=M.DSTable(8)
def pred(W,w):
    T=T8
    return sum(F(c)*monval(m) for c,m in zip(T.reduce_word(w),T.monomials))
# depth-2 weight-8 words vs dz at >=60 digits
print('--- depth-2 W8 vs dz ---')
for w in [(2,6),(3,5),(4,4),(5,3),(6,2),(7,1)]:
    a,b=w; val=dz(a,b); pr=pred(8,w)
    print(w, 'agree digits:', nstr(abs(val-pr),3))
print('t=%.1fs'%(time.time()-t0))
