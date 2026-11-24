import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=80
a1v=1e7
t0=time.time()
# detailed scan: track sheets and record Q11 sign per sheet
prev=None; changes=[]
a2grid=np.logspace(-4,9,140)
for a2v in a2grid:
    cur=[]
    for (a4v,a5v) in v10_mpf(a1v,a2v):
        cur.append((a4v,a5v,q11_mpf([a1v,a2v,1,a4v,a5v])))
    if prev:
        a2p,prevl=prev
        for (a4v,a5v,q) in cur:
            if not prevl: continue
            j=min(range(len(prevl)), key=lambda k:(prevl[k][0]-a4v)**2/max(1,a4v**2)+(prevl[k][1]-a5v)**2/max(1,a5v**2))
            a4p,a5p,qp=prevl[j]
            if qp*q<0:
                changes.append((a2p,a2v,a4v,a5v,qp,q))
    prev=(a2v,cur)
print(f"a1=1e7: {len(changes)} sign-change brackets [{time.time()-t0:.0f}s]")
for (a2p,a2v,a4v,a5v,qp,q) in changes[:20]:
    print(f"  a2 in [{a2p:.4g},{a2v:.4g}], a4={mp.nstr(a4v,6)} a5={mp.nstr(a5v,6)}, Q11 {mp.nstr(qp,2)}->{mp.nstr(q,2)}")
