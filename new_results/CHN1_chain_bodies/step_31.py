import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import engine as E
# engine confirmation of plane vanishing
for a in [(2,2,3,3),(1,1,2,2),(2,3,2,3),(3,1,3,1)]:
    r=E.reduced_coeffs(a,8)[8]
    print(a, 'z(3,5) coeff =', r[3], 'plane', 'Pi1' if a[0]==a[1] and a[2]==a[3] else ('Pi2' if a[0]==a[2] and a[1]==a[3] else '?'))
# blind direction prediction
a=(3,1,1,1)
r8=E.reduced_coeffs(a,8)[8]
print('blind direction', a, 'c8 coords:', r8)
open('progress.log','a').write('step_26-27: Q4 = (a1+a2+a3+a4)*C3/10, C3 irreducible / smooth cubic; vanishes on Pi1,Pi2 (engine-confirmed); grad C3(1,1,1,1)=27(1,-1,-1,1); blind dir (3,1,1,1) chosen\n')
print('t=%.1fs'%(time.time()-t0))
