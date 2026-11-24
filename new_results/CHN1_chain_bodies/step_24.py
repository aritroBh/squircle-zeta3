import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
open('progress.log','a').write('step_21: control C PASS: exact c2..c8 vs all 6 pilot files; diag diffs 5e-51..3e-49 (<1e-45), generic 1e-70..4.6e-65 (<1e-60)\n')
t0=time.time()
from mpmath import mp, mpf, quad, gamma, hyp3f2, nstr, pi
mp.dps=35
def xintegral(us):
    u1,u2,u3,u4=us; p1,p2,p3,p4=[1/x for x in us]
    return quad(lambda a: (1-a**p2)**u1*quad(lambda b: (1-b**p3)**u4, [0,(1-a**p2)**(1/p3)]), [0,mpf(1)/2,mpf('0.9'),1])
def formula(us):
    u1,u2,u3,u4=us
    return gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
pts=[(mpf(1)/2,mpf(1)/3,mpf(1)/4,mpf(1)/5),
     (mpf(1),mpf(1)/2,mpf(1)/3,mpf(2)/3),
     (mpf(2)/3,mpf(1)/2,mpf(2)/5,mpf(1)/2)]
for us in pts:
    t=time.time()
    d=xintegral(us); F=formula(us)
    print([str(x) for x in us], '|d-F|=',nstr(abs(d-F),3), '(%.0fs)'%(time.time()-t))
# Vol(1/2)=4+pi^2/2
us=(mpf(1)/2,)*4
t=time.time()
F=formula(us)
print('Vol(1/2)*16 vs 4+pi^2/2: |diff|=', nstr(abs(16*F-(4+pi**2/2)),3), '(%.0fs)'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
