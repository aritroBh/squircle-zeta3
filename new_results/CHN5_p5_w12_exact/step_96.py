import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=70
import math
def dense_count(a1v):
    """independent: dense a2 scan (wide range), sign-change -> Newton; count + signs. No continuation."""
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    lo=math.log10(max(a1v,1e-7))-4; hi=math.log10(max(a1v,1e-7))+3.5
    prev=None; guesses=[]
    for a2v in np.logspace(lo,hi,90):
        cur=[]
        for (a4v,a5v) in v10_mpf(a1v,a2v):
            cur.append((a4v,a5v,q11_mpf([a1v,a2v,1,a4v,a5v])))
        if prev:
            a2p,prevl=prev
            for (a4v,a5v,q) in cur:
                if not prevl: continue
                j=min(range(len(prevl)), key=lambda k:(prevl[k][0]-a4v)**2/max(1,a4v**2)+(prevl[k][1]-a5v)**2/max(1,a5v**2))
                a4p,a5p,qp=prevl[j]
                if qp*q<0: guesses.append(((a2p+a2v)/2,(a4p+a4v)/2,(a5p+a5v)/2))
        prev=(a2v,cur)
    sols=[]
    for g in guesses:
        try:
            r=mp.findroot(f, tuple(mp.mpf(x) for x in g), tol=mp.mpf('1e-40'), maxsteps=200)
            a2r,a4r,a5r=r[0],r[1],r[2]
            if a2r>0 and a4r>0 and a5r>mp.mpf('1e-25'):
                key=tuple(mp.nstr(x,15) for x in (a2r,a4r,a5r))
                if not any(key==k for k,_ in sols): sols.append((key,(a2r,a4r,a5r)))
        except Exception: pass
    out=[]
    for _,(a2r,a4r,a5r) in sols:
        pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
        out.append(('+' if mpoly(QA,pt)>0 else '-', '+' if mpoly(QC,pt)>0 else '-'))
    return out
t0=time.time()
for a1v in [1e-5, 0.003, 0.05, 0.4, 0.8, 1.5, 20, 300, 5000, 2e5, 3e6]:
    out=dense_count(a1v)
    print(f"a1={a1v:.3g}: {len(out)} interior Gamma pts, signs(QA,QC)={out} [{time.time()-t0:.0f}s]")
