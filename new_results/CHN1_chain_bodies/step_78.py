import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
open('progress.log','a').write('step_70-74: products span dim 10 at W12; diag 5c10 matches pilot; C3xQ5 = 14 lines, positive common zeros EXACTLY Pi1 u Pi2; Q5 grad at diag = 2988(1,-1,-1,1); pilot c9,c10,c11 all match within pilot accuracy\n')
src=r'''
"""Diagonal c12 Cauchy validation: a=(1,1,1,1), R=1/3, r=1/18, N=128, dps 85.
Predicted: c12 from exact engine (products + z39/z57 only; w*=0)."""
import os, json, time, sys
from mpmath import mp, mpf, mpc, gamma, quad, hyp2f1, expjpi, nstr
os.chdir('<project>/research/08_theory/CHN1_chain_bodies')
mp.dps=85
A=(1,1,1,1); N=128; R=mpf(1)/3; r=R/6
def Vol16(u):
    u1,u2,u3,u4=[a*u for a in A]
    F3 = 1 + u3*quad(lambda tt: tt**(u3-1)*(hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt)-1), [0,1])
    return gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*F3
pts={}
for j in range(0, N//2+1):
    u=r*expjpi(mpf(2*j)/N)
    pts[j]=Vol16(u)
vals=[pts[j] if j<=N//2 else pts[N-j].conjugate() for j in range(N)]
res={}
for k in range(10,15):
    ck=sum(vals[j]*expjpi(mpf(-2*j*k)/N) for j in range(N))/N/r**k
    res[k]=(nstr(ck.real,75), nstr(ck.imag,3))
json.dump(res, open('diag_c12_result.json','w'), indent=1)
print('c12 =', res[12])
'''
open('diag_c12.py','w').write(src)
import subprocess, sys
proc=subprocess.Popen(['nice','-n','10','nohup',sys.executable,'diag_c12.py','>','diag_c12.out','2>&1'],
                      cwd=WORK, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print('launched', proc.pid)
t0=time.time()
import mzv_exact as M, engine as E
from fractions import Fraction as F
def table(W):
    if W not in E._TABLES: E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
labels=['z12','z3z9','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z35','z2z37','z39','z57','w*']
for a in [(1,1,1,2),(1,2,2,1)]:
    r=E.reduced_coeffs(a,12)[12]
    print(a, {labels[i]:r[i] for i in (9,10,11,7,8)})
print('t=%.1fs'%(time.time()-t0))
