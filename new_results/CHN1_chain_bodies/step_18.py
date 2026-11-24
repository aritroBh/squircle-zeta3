import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
src=r'''
"""Exact engine: Taylor coefficients c_W(a) of Vol/16 along u_i = a_i u, as
Q-combinations of admissible MZV words of weight W. Then reduced to the DS basis.
Vol/16 = G * 3F2(-u4,u3,1+u1+u3; 1+u3,1+u1+u2+u3; 1), u_i = a_i u."""
from fractions import Fraction as F
import mzv_exact as M

def st_letter(i, d):
    """stuffle single letter (i,) into every word of dict d."""
    nd={}
    for w,c in d.items():
        for w2,c2 in M.stuffle((i,),w).items():
            nd[w2]=nd.get(w2,F(0))+c*c2
    return {w:c for w,c in nd.items() if c}

def st_dicts(d1,d2):
    nd={}
    for w1,c1 in d1.items():
        for w2,c2 in d2.items():
            for w,c in M.stuffle(w1,w2).items():
                nd[w]=nd.get(w,F(0))+c1*c2*c
    return {w:c for w,c in nd.items() if c}

def coeffs(a, Wmax=11):
    a1,a2,a3,a4=a
    A=a1+a3; B=a1+a2+a3
    d=[F(0)]*(Wmax+2)
    g=[F(0)]*(Wmax+2)
    L=[F(0)]*(Wmax+2)
    for m in range(1,Wmax+2):
        sgn = 1 if m%2==1 else -1   # (-1)^(m+1)
        d[m]=F(sgn*(A**m-B**m), m)
        g[m]=d[m]-F(a4**m, m)
        L[m]=F((-1 if m%2==0 else 1)*(a2**m + A**m - B**m), m)  # (-1)^m/m (...)
    # E3: delta
    delta=[F(0)]*(Wmax+1); delta[0]=F(1)
    for l in range(1,Wmax+1):
        delta[l]=F(1,l)*sum(F(i)*d[i]*delta[l-i] for i in range(1,l+1))
    # F[j]: exp(sum g_m u^m S_m(k-1)) -> word dicts
    Fw=[{} for _ in range(Wmax+1)]; Fw[0]={():F(1)}
    for j in range(1,Wmax+1):
        acc={}
        for i in range(1,j+1):
            acc=M.wd_add_scaled(acc, st_letter(i, Fw[j-i]), F(i*g[i],j))
        Fw[j]=acc
    # H
    H=[{} for _ in range(Wmax+1)]; H[0]={():F(1)}
    for W in range(2,Wmax+1):
        acc={}
        for r in range(0,W-1):
            for l in range(0,W-1-r):
                j=W-2-r-l
                c0=-a3*a4*((-a3)**r)*delta[l]
                if c0==0: continue
                s=2+r+l
                for w,c in Fw[j].items():
                    acc[(s,)+w]=acc.get((s,)+w,F(0))+c0*c
        H[W]={w:c for w,c in acc.items() if c}
    # G
    G=[{} for _ in range(Wmax+1)]; G[0]={():F(1)}
    for W in range(2,Wmax+1):
        acc={}
        for i in range(2,W+1):
            acc=M.wd_add_scaled(acc, st_letter(i, G[W-i]), F(i*L[i],W))
        G[W]=acc
    # c = G*H
    c=[{} for _ in range(Wmax+1)]
    c[0]={():F(1)}
    for W in range(1,Wmax+1):
        acc={}
        for i in range(0,W+1):
            acc=M.wd_add_scaled(acc, st_dicts(G[i],H[W-i]), F(1))
        c[W]={w:cc for w,cc in acc.items() if cc}
    # assertions: weight and admissibility
    for W in range(1,Wmax+1):
        for w in c[W]:
            assert sum(w)==W and w[0]>=2, (W,w)
    return c

_TABLES={}
def table(W):
    if W not in _TABLES: _TABLES[W]=M.DSTable(W)
    return _TABLES[W]

def reduced_coeffs(a, Wmax=11):
    c=coeffs(a,Wmax)
    out={}
    for W in range(2,Wmax+1):
        T=table(W)
        acc=[F(0)]*T.nullity
        for w,cc in c[W].items():
            rw=T.reduce_word(w)
            for i in range(T.nullity): acc[i]+=cc*rw[i]
        out[W]=acc
    return out
'''
open('engine.py','w').write(src)
print('engine.py written')
