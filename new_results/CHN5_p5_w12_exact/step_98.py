import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=80
t0=time.time()
# small-a1 asymptotic: solve at a1=1e-8..1e-11 near the known branch (a2~2.6667, a4~2.6667, a5~1)
for a1v in [1e-8, 1e-9, 1e-10, 1e-11]:
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    r=mp.findroot(f, (mp.mpf('2.6667'),mp.mpf('2.6667'),mp.mpf('1.0')), tol=mp.mpf('1e-55'), maxsteps=300)
    a2r,a4r,a5r=r[0],r[1],r[2]
    pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
    S=sum(pt)
    qa=mpoly(QA,pt); qc=mpoly(QC,pt)
    print(f"a1={a1v:.0e}: a2={mp.nstr(a2r,10)} a4={mp.nstr(a4r,10)} a5={mp.nstr(a5r,10)}")
    print(f"   QA/S^8={mp.nstr(qa/S**8,6)} QC/S^8={mp.nstr(qc/S**8,6)}  log10|QA|/log10(a1)={mp.log10(abs(qa/S**8))/mp.log10(a1v if a1v<1 else 1):.4f}")
print(f"[{time.time()-t0:.0f}s]")
