import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, nstr, gamma, hyp3f2, quad, hyp2f1
mp.dps=85
u=mpf(1)/60
# sanity 1: integral of t^{u-1}
v=quad(lambda tt: tt**(u-1), [0,1])
print('int t^{u-1}: err', nstr(abs(v-1/u),3))
v=quad(lambda tt: tt**(u-1), [0,1], maxdegree=10)
print('int t^{u-1} maxdegree=10: err', nstr(abs(v-1/u),3))
# sanity 2: the 3F2 identity at u=1/2
u2=mpf(1)/2
F3 = u2*quad(lambda tt: tt**(u2-1)*hyp2f1(-u2,1+4*u2,1+5*u2,tt), [0,1])
H = hyp3f2(-u2,u2,1+4*u2,1+u2,1+5*u2,1)
print('identity u=1/2: |diff|', nstr(abs(F3-H),3))
print('t=%.1fs'%(time.time()-t0))
