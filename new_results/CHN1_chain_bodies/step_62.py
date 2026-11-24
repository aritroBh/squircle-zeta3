import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
from fractions import Fraction as F
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
# sanity: identity and dependent rows
M1=[[F(1),F(0)],[F(0),F(1)],[F(1),F(1)]]
print('rank test (expect 2):', rank(M1))
M2=[[F(1),F(2)],[F(2),F(4)],[F(1),F(1)]]
print('rank test (expect 2):', rank(M2))
# 11 independent + 1 dependent
import itertools
M3=[[F(1) if i==j else F(0) for j in range(12)] for i in range(11)]+[[F(0)]*12]
print('rank test (expect 11):', rank(M3))
