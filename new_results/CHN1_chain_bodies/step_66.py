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
wstar=(2,1,1,8)
print('dual of (2,1,1,1,1,1,1,4):', M.dual((2,1,1,1,1,1,1,4)))
print('dual of wstar:', M.dual(wstar))
# duality at W=12 (table check)
dw=all(elim({w:F(1)})==elim({M.dual(w):F(1)}) for w in M.WORDS[12])
print('duality all W=12 words:', dw)
# dependency relation among old 12 monomials
red=[vec(elim(M.wordcombo(m))) for m in M.BASIS[12]]
# find kernel vector: solve sum c_i red[i]=0
A=[r[:] for r in red]
# rref on transposed: find combination
import itertools
Mtx=[[red[j][i] for j in range(12)] for i in range(12)]
# solve Mtx c = 0
n=12
B=[row[:] for row in Mtx]
pivcols=[];r0=0
for c in range(n):
    p=next((r for r in range(r0,n) if B[r][c]!=0),None)
    if p is None: continue
    B[r0],B[p]=B[p],B[r0]
    inv=F(1)/B[r0][c]; B[r0]=[x*inv for x in B[r0]]
    for r in range(n):
        if r!=r0 and B[r][c]!=0:
            f=B[r][c]; B[r]=[x-f*y for x,y in zip(B[r],B[r0])]
    pivcols.append(c); r0+=1
freec=[c for c in range(n) if c not in pivcols]
labels=['z12','z3z9','z5z7','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z35','z39','z57','z2z37']
for fc in freec:
    cvec=[F(0)]*n; cvec[fc]=F(1)
    for k,c in enumerate(pivcols):
        cvec[c]=-B[k][fc] if k<len(B) else F(0)
    print('dependency:', [(labels[i],cvec[i]) for i in range(n) if cvec[i]!=0])
# candidate basis: 7 products + z12,z3z9,z5z7,z39 + wstar ; check rank 12
prods_idx=[3,4,5,6,7,8,11]
cand=[red[i] for i in prods_idx]+[red[0],red[1],red[2],red[9]]+[vec(elim({wstar:F(1)}))]
print('candidate basis rank:', rank(cand))
cand2=[red[i] for i in prods_idx]+[red[0],red[1],red[2],red[10]]+[vec(elim({wstar:F(1)}))]
print('candidate basis (with z57 instead of z39) rank:', rank(cand2))
print('t=%.1fs'%(time.time()-t0))
