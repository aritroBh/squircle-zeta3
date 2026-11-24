import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
# numeric validation of the W=12 relation first
from mpmath import mp, mpf, zeta, nstr
mp.dps=60
from mzv_vals import dz
lhs = dz(5,7)
rhs = mpf(128195)/12438*zeta(12) - mpf(14)/9*zeta(3)*zeta(9) - mpf(25)/3*zeta(5)*zeta(7) + mpf(14)/9*dz(3,9)
print('W12 relation residual:', nstr(abs(lhs-rhs),3))
# also numeric check z(3,9) reduction exists in span: not needed further
# update BASIS[12]
import mzv_exact as M
newb=[((12,),),((3,),(9,)),((2,),(3,),(7,)),((2,),(5,),(5,)),((2,),(2,),(3,),(5,)),
      ((2,),(2,),(2,),(3,),(3,)),((3,),(3,),(3,),(3,)),((2,),(2,),(3,5)),((2,),(3,7)),
      ((3,9),),((5,7),),((2,1,1,8),)]
src=open('mzv_exact.py').read()
old="""BASIS[12]=[((12,),),((3,),(9,)),((5,),(7,)),((2,),(3,),(7,)),((2,),(5,),(5,)),((2,),(2,),(3,),(5,)),
           ((2,),(2,),(2,),(3,),(3,)),((3,),(3,),(3,),(3,)),((2,),(2,),(3,5)),((3,9),),((5,7),),((2,),(3,7))]"""
new="""BASIS[12]=[((12,),),((3,),(9,)),((2,),(3,),(7,)),((2,),(5,),(5,)),((2,),(2,),(3,),(5,)),
           ((2,),(2,),(2,),(3,),(3,)),((3,),(3,),(3,),(3,)),((2,),(2,),(3,5)),((2,),(3,7)),
           ((3,9),),((5,7),),((2,1,1,8),)]"""
assert old in src
src=src.replace(old,new)
open('mzv_exact.py','w').write(src)
print('BASIS[12] updated')
print('t=%.1fs'%(time.time()-t0))
