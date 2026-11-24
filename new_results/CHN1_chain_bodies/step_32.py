import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr, gamma, hyp3f2, quad, hyp2f1, mpc, expjpi
mp.dps=85
from mzv_vals import dz
from fractions import Fraction as F
import engine as E
z35=dz(3,5,N=800,J=60)
r8=[F(-1421131,40),F(104724,5),F(-2970),F(1746,5)]
pred = (mpf(int(r8[0].numerator))/r8[0].denominator)*zeta(8) + (mpf(int(r8[1].numerator))/r8[1].denominator)*zeta(3)*zeta(5) + (mpf(int(r8[2].numerator))/r8[2].denominator)*zeta(2)*zeta(3)**2 + (mpf(int(r8[3].numerator))/r8[3].denominator)*z35
open('blind_prediction.txt','w').write(nstr(pred,80)+'\n')
print('PREDICTED c8 at (3,1,1,1) =', nstr(pred,50))
# time hyp3f2 at a complex u
u = mpc(mpf(1)/30)*expjpi(mpf(2)/7)
t=time.time()
u1,u2,u3,u4=3*u,u,u,u
v1 = gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
print('hyp3f2 path: %.0fs'%(time.time()-t))
# integral path
t=time.time()
F3 = u3*quad(lambda tt: tt**(u3-1)*hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt), [0,1])
v2 = gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*F3
print('integral path: %.0fs'%(time.time()-t))
print('|v1-v2| =', nstr(abs(v1-v2),3))
print('t=%.1fs'%(time.time()-t0))
