import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from mpmath import mp, mpf, zeta, nstr, gamma, hyp3f2, quad, hyp2f1, mpc
mp.dps=85
u=mpf(1)/60
u1,u2,u3,u4=3*u,u,u,u
t=time.time()
v1 = gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
print('hyp3f2 real: %.0fs'%(time.time()-t))
t=time.time()
F3 = u3*quad(lambda tt: tt**(u3-1)*hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt), [0,1])
v2 = gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*F3
print('quad default: %.0fs  |v1-v2|=%s'%(time.time()-t, nstr(abs(v1-v2),3)))
# manual tanh-sinh with more levels
from mpmath.calculus.quadrature import QuadratureRule, TanhSinh
ts=TanhSinh(mp)
t=time.time()
F3b = u3*ts.sumfun(lambda tt: tt**(u3-1)*hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt), 0, 1, tol=mp.eps*100, maxnodes=10**6)
print('manual ts: %.0fs'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
