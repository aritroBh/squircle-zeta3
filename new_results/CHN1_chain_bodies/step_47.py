import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import itertools
# hardcode the Q_m printed (verified from previous step's output)
def solve_forms(Q, cands, mrange):
    # Q: dict m -> dict (i,j)->Fraction coefficient of s^i t^j ; cands: list of (a,b)
    # equations: for each m, each monomial: sum n_c * monomial_coeff((a s + b t)^m) = Q_m coeff
    from math import comb
    rows=[]; rhs=[]
    for m in mrange:
        for i in range(m+1):
            j=m-i
            row=[F(comb(m,i))*F(a)**i*F(b)**j for (a,b) in cands]
            rows.append(row)
            rhs.append(Q.get(m,{}).get((i,j),F(0)))
    # exact gaussian elimination, find particular solution
    nr=len(rows); nc=len(cands)
    A=[r[:]+[b] for r,b in zip(rows,rhs)]
    piv=[]; r0=0
    for c in range(nc):
        p=None
        for r in range(r0,nr):
            if A[r][c]!=0: p=r; break
        if p is None: continue
        A[r0],A[p]=A[p],A[r0]
        inv=F(1)/A[r0][c]; A[r0]=[x*inv for x in A[r0]]
        for r in range(nr):
            if r!=r0 and A[r][c]!=0:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[r0])]
        piv.append(c); r0+=1
        if r0==nr: break
    # consistency
    for r in range(r0,nr):
        if any(A[r][c]!=0 for c in range(nc)): continue
        if A[r][nc]!=0: return None,'INCONSISTENT'
    free=[c for c in range(nc) if c not in piv]
    sol=[F(0)]*nc
    for k,c in enumerate(piv): sol[c]=A[k][nc]
    return (sol, free, piv, A), 'OK'
# test on Pi2 Q data
def Q_poly(entries):
    # entries: dict m -> list of (i, coeff) for s^i t^(m-i)
    out={}
    for m,lst in entries.items():
        out[m]={}
        for i,c in lst: out[m][(i,m-i)]=F(c)
    return out
Q_Pi2 = Q_poly({
 2:[(1,-6)],
 3:[(1,-15),(2,-15)],
 4:[(1,-36),(2,-4),(3,-36)],
 5:[(1,-85),(2,-170),(3,-170),(4,-85)],
 6:[(1,-198),(2,-138),(3,-282),(4,-138),(5,-198)],
 7:[(1,-455),(2,-1365),(3,-2275),(4,-2730),(5,-2275),(6,-1365),(7,-455)],
 8:[(1,-1032),(2,-1590),(3,-3668),(4,-3916),(5,-3916),(6,-3668),(7,-1590),(8,-1032)],   # -8st(387s6+477s5t+979s4t2+995s3t3+979s2t4+477st5+387t6)/3: coeffs: -8*387/3=-1032, -8*477/3=-1272? check below
})
print('placeholder - need exact Q8..Q11; recompute directly instead')
print('t=%.1fs'%(time.time()-t0))
