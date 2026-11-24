import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M
from fractions import Fraction as F
rels=M.relations(12)
idx={w:i for i,w in enumerate(M.WORDS[12])}
piv=M.rref_flint(rels, idx)
free=[j for j in range(len(idx)) if all(c!=j for c,_ in piv)]
def elim(v):
    r={idx[w]:F(c) for w,c in v.items()}
    for col,pr in piv:
        c=r.get(col)
        if c:
            for j,pc in pr.items():
                r[j]=r.get(j,F(0))-c*pc
                if r[j]==0: del r[j]
    return r
def vec(v): return [v.get(fc,F(0)) for fc in free]
def rank(rows):
    A=[list(r) for r in rows]; r0=0; n=len(A[0])
    for c in range(n):
        p=next((r for r in range(r0,len(A)) if A[r][c]!=0),None)
        if p is None: continue
        A[r0],A[p]=A[p],A[r0]
        inv=F(1)/A[r0][c]; A[r0]=[x*inv for x in A[r0]]
        for r in range(len(A)):
            if r!=r0 and A[r][c]!=0:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[r0])]
        r0+=1
    return r0
red=[vec(elim(M.wordcombo(m))) for m in M.BASIS[12]]
rk=rank(red)
print('span of 12 monomials:', rk)
found=None
for w in M.WORDS[12]:
    if len(w)<3: continue
    v=vec(elim({w:F(1)}))
    if rank(red+[v])>rk:
        found=w; break
print('a word outside the span:', found)
print('t=%.1fs'%(time.time()-t0))
