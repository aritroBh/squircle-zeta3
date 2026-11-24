import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import pickle, sympy as sp
from fractions import Fraction as F
from math import comb
cache=pickle.load(open('locus_logV.pkl','rb'))
s,t=sp.symbols('s t')
def Q_to_dicts(name):
    out={}
    for m,poly in cache[name]['Q'].items():
        P=sp.Poly(poly, s, t)
        out[m]={}
        for mon,c in P.terms(): out[m][mon]=F(c.numerator, c.denominator)
    return out
def solve_forms(Q, cands, mrange, rhs_override=None):
    rows=[]; rhs=[]
    for m in mrange:
        for i in range(m+1):
            j=m-i
            rows.append([F(comb(m,i))*(F(a)**i)*(F(b)**j) for (a,b) in cands])
            rhs.append(Q.get(m,{}).get((i,j),F(0)))
    nr=len(rows); nc=len(cands)
    A=[r[:]+[b] for r,b in zip(rows,rhs)]
    piv=[]; r0=0
    for c in range(nc):
        p=next((r for r in range(r0,nr) if A[r][c]!=0), None)
        if p is None: continue
        A[r0],A[p]=A[p],A[r0]
        inv=F(1)/A[r0][c]; A[r0]=[xx*inv for xx in A[r0]]
        for r in range(nr):
            if r!=r0 and A[r][c]!=0:
                f=A[r][c]; A[r]=[xx-f*yy for xx,yy in zip(A[r],A[r0])]
        piv.append(c); r0+=1
        if r0==nr: break
    for r in range(r0,nr):
        if A[r][nc]!=0 and all(A[r][c]==0 for c in range(nc)):
            return None, None
    free=[c for c in range(nc) if c not in piv]
    sol=[F(0)]*nc
    for k,c in enumerate(piv): sol[c]=A[k][nc]
    return (sol,free,piv), (rows,rhs)
for name in ['Pi2','Pi1']:
    Q=Q_to_dicts(name)
    res,_=solve_forms(Q, cands:=[(a,b) for a in range(-3,4) for b in range(-3,4) if (a,b)!=(0,0)], range(2,9))
    print('=== %s (m<=8) ==='%(name))
    if res:
        s0,free,piv=res
        nz=[(cands[i],s0[i]) for i in range(len(cands)) if s0[i]!=0]
        print('  sol:', nz, ' free:', [cands[i] for i in free])
print('t=%.1fs'%(time.time()-t0))
