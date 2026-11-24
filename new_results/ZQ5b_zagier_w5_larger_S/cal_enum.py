
import sys, time
import numpy as np
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F = W.make_field(13, split=[3], inert=[2])
gens=F['letters']; d=len(gens)
Lf = np.array([float(mp.log(abs(g.emb(60))) if False else mp.log(g.emb(60)) if g.emb(60)>0 else mp.log(-g.emb(60))) for g in gens])
B=8
axes=[np.arange(-B,B+1)]*d
grids=np.meshgrid(*axes,indexing='ij')
T=np.stack([g.ravel() for g in grids],axis=1).astype(float)
logabs = T @ Lf
mask = logabs < 1e-9
keep = T[mask].astype(int)
keep = keep[np.any(keep!=0,axis=1)]
print(f"d={d} B={B}: {len(T)} tuples, {len(keep)} with |x|<1")
# time exact build + factor on first 20000
t0=time.time()
found=0
for exps in keep[:20000]:
    xp=QD(1,0,13)
    for g,e in zip(gens,exps):
        if e: xp=xp*(g**int(e))
    v=xp.emb(70)
    if not (abs(v)<1): continue
    for sgn in (1,-1):
        x=xp if sgn==1 else -xp
        f=W.factor_sunit(1-x,F)
        if f is not None: found+=1
dt=time.time()-t0
print(f"20000 tuples: {dt:.1f}s, found {found} -> full run ~ {dt*len(keep)/20000:.0f}s")
