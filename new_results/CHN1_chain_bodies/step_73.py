import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
open('progress.log','a').write('step_66-67: W12 PLANES PROVEN: obs1=z39+14/9 z57 = 0 and w*=0,z2^2z35=0,z2z37=0 at 15 pts/plane (exact) => c12 in products span on Pi1,Pi2,diag. Thomae orbit(120) on Pi1,Pi2 has Saalschutz elements; diag orbit NONE\n')
t0=time.time()
from fractions import Fraction as F
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES: E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
# products span rank at W=12
prods=[((12,),),((3,),(9,)),((5,),(7,)),((2,),(3,),(7,)),((2,),(5,),(5,)),((2,),(2,),(3,),(5,)),
      ((2,),(2,),(2,),(3,),(3,)),((3,),(3,),(3,),(3,)),((2,),(2,),(3,5)),((2,),(3,7))]
T12=table(12)
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
# value vectors of products in quotient: reduce each onto free cols
def elimvec(v):
    r=T12._elim(v); return [r.get(fc,F(0)) for fc in T12.freecols]
rv=[elimvec(M.wordcombo(m)) for m in prods]
print('W=12 products span rank:', rank(rv), '(of 10; quotient dim = 12-rank)')
# diagonal c10 vs pilot vector
r10=E.reduced_coeffs((1,1,1,1),10)[10]
print('diag 5c10 =', [5*x for x in r10], ' expect [-85307, 92520, 46260, -30960, -11560, 0, 0]')
print('t=%.1fs'%(time.time()-t0))
