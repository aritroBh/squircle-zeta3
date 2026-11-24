import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=80
t0=time.time()
# check Q11 signs on all V10 sheets at a1=1e7 and a1=3e6
for a1v in [3e6, 1e7]:
    q11s=[]
    for a2v in np.logspace(-4,8,100):
        for (a4v,a5v) in v10_mpf(a1v,a2v):
            q=q11_mpf([a1v,a2v,1,a4v,a5v])
            q11s.append(q)
    sgns=set('+' if q>0 else '-' for q in q11s)
    mn=min((abs(q) for q in q11s), default=None)
    print(f"a1={a1v:.3g}: {len(q11s)} V10 sheet pts, Q11 signs={sgns}, min|Q11|={mp.nstr(mn,4) if mn else 'NA'} [{time.time()-t0:.0f}s]")
