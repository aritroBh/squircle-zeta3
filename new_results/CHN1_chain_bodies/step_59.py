import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M
from fractions import Fraction as F
rels=M.relations(12)
idx={w:i for i,w in enumerate(M.WORDS[12])}
t=time.time()
piv=M.rref_flint(rels, idx)
print('rank(relations)=%d words=%d nullity=%d (%.0fs)'%(len(piv),len(idx),len(idx)-len(piv),time.time()-t))
free=[j for j in range(len(idx)) if all(c!=j for c,_ in piv)]
# reduce each monomial onto free cols
def elim(v):
    r={idx[w]:F(c) for w,c in v.items()}
    for col,pr in piv:
        c=r.get(col)
        if c:
            for j,pc in pr.items():
                r[j]=r.get(j,F(0))-c*pc
                if r[j]==0: del r[j]
    return r
mv=[M.wordcombo(m) for m in M.BASIS[12]]
red=[elim(v) for v in mv]
# rank of reduced monomials
rows=[[rr.get(fc,F(0)) for fc in free] for rr in red]
import itertools
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
    return r0, A
rk,A=rank(rows)
print('monomial rank:', rk, 'of', len(rows))
labels=['z12','z3z9','z5z7','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z35','z39','z57','z2z37']
# find dependency: which row becomes zero
for i,r in enumerate(A):
    if all(x==0 for x in r): print('dependent:', labels[i])
print('t=%.1fs'%(time.time()-t0))
