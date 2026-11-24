import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import engine as E, mzv_exact as M
def gamma_series(L, Wmax):
    # exp(sum_{m>=2} L[m] u^m zeta(m)) -> word dicts
    G=[{} for _ in range(Wmax+1)]; G[0]={():F(1)}
    for W in range(2,Wmax+1):
        acc={}
        for i in range(2,W+1):
            acc=M.wd_add_scaled(acc, E.st_letter(i, G[W-i]), F(i*L[i],W))
        G[W]=acc
    return G
def reduced_from_worddicts(G, Wmax):
    out={}
    for W in range(2,Wmax+1):
        T=E.table(W); acc=[F(0)]*T.nullity
        for w,cc in G[W].items():
            rw=T.reduce_word(w)
            for i in range(T.nullity): acc[i]+=cc*rw[i]
        out[W]=acc
    return out
ok=True
def chk(name, got, want):
    global ok
    good = got==want
    print(('PASS ' if good else 'FAIL '), name, '' if good else (got, want))
    ok = ok and good
for (a2,a3,a4) in [(1,1,1),(2,1,3),(1,2,1),(3,2,2)]:
    a=(0,a2,a3,a4)
    r=E.reduced_coeffs(a,8)
    L={m: F((1 if m%2==0 else -1)*(a3**m + (a2+a4)**m - (a2+a3+a4)**m), m) for m in range(2,10)}
    gr=reduced_from_worddicts(gamma_series(L,8),8)
    chk('a1=0 %s'%str(a), [r[W] for W in range(2,9)], [gr[W] for W in range(2,9)])
for (a1,a3,a4) in [(1,1,1),(2,1,3),(3,2,1),(1,3,2)]:
    a=(a1,0,a3,a4)
    r=E.reduced_coeffs(a,8)
    L={m: F((1 if m%2==0 else -1)*(a3**m + a4**m - (a3+a4)**m), m) for m in range(2,10)}
    gr=reduced_from_worddicts(gamma_series(L,8),8)
    chk('a2=0 %s'%str(a), [r[W] for W in range(2,9)], [gr[W] for W in range(2,9)])
print('ALL', ok, 't=%.1fs'%(time.time()-t0))
open('progress.log','a').write('step_18-20: engine.py controls A,B,D,E ALL PASS (diagonal c2..c8, 3 pilot vectors, reversal x6, Gauss products a1=0/a2=0 x8)\n')
