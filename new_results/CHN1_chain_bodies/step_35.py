import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, mpc, nstr, gamma, hyp3f2, quad, hyp2f1, expjpi
mp.dps=85
# diagnose the plain power integral
mp.dps=50
v=quad(lambda tt: tt**(mpf(1)/60-1), [0,mpf('0.001'),1])
print('split int t^(u-1): err', nstr(abs(v-60),3))
mp.dps=85
def F3_sub(u1,u2,u3,u4):
    return 1 + u3*quad(lambda tt: tt**(u3-1)*(hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt)-1), [0,1])
# real u
u=mpf(1)/60; u1,u2,u3,u4=3*u,u,u,u
v1 = hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
t=time.time(); v2=F3_sub(u1,u2,u3,u4)
print('real u: |sub-hyp3f2| =', nstr(abs(v1-v2),3), '(%.0fs)'%(time.time()-t))
# complex u with Re u3 < 0 (the continuation regime)
uc = mpc(mpf(1)/30)*expjpi(mpf(2)*3/7)
u1,u2,u3,u4=3*uc,uc,uc,uc
print('Re u3 =', nstr(uc.real,4))
t=time.time()
v1 = hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
print('hyp3f2 complex: %.0fs'%(time.time()-t))
t=time.time()
v2=F3_sub(u1,u2,u3,u4)
print('complex u: |sub-hyp3f2| =', nstr(abs(v1-v2),3), '(%.0fs)'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
