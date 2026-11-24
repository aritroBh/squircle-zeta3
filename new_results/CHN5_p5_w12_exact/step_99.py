import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=80
t0=time.time()
for a1v in [1e-7, 1e-8, 1e-9, 1e-10, 1e-11]:
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    r=mp.findroot(f, (mp.mpf('2.6667'),mp.mpf('2.6667'),mp.mpf('1.0')), tol=mp.mpf('1e-55'), maxsteps=300)
    a2r,a4r,a5r=r[0],r[1],r[2]
    pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
    S=sum(pt)
    qa=mpoly(QA,pt); qc=mpoly(QC,pt)
    nqa=qa/S**8; nqc=qc/S**8
    print(f"a1={a1v:.0e}: a2={mp.nstr(a2r,9)} a5={mp.nstr(a5r,9)}  nQA={mp.nstr(nqa,6)} nQC={mp.nstr(nqc,6)}  nQA/a1={mp.nstr(nqa/a1v,6)}")
print(f"[{time.time()-t0:.0f}s]")
