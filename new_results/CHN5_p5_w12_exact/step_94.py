import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=80
a1v=1e7
C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
brackets=[(6e6,9.5e6,3.46e6),(1.87e7,2.07e7,2.74e6),(2.88e7,2.93e7,3.2e7),(6.8e7,7.57e7,3.65e7)]
t0=time.time()
for g in brackets:
    try:
        r=mp.findroot(f, tuple(mp.mpf(x) for x in g), tol=mp.mpf('1e-45'), maxsteps=300)
        a2r,a4r,a5r=r[0],r[1],r[2]
        pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
        S=sum(pt)
        qa=mpoly(QA,pt); qc=mpoly(QC,pt)
        res=max(abs(fi(a2r,a4r,a5r)) for fi in f)
        print(f"a2={mp.nstr(a2r,10)} a4={mp.nstr(a4r,10)} a5={mp.nstr(a5r,10)}")
        print(f"   QA={'+' if qa>0 else '-'} QC={'+' if qc>0 else '-'}  QA/S^8={mp.nstr(qa/S**8,4)} QC/S^8={mp.nstr(qc/S**8,4)}  res={mp.nstr(res,2)}")
    except Exception as e:
        print("fail", e)
print(f"[{time.time()-t0:.0f}s]")
