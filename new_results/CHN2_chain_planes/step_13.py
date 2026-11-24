
import pickle, time, json
from fractions import Fraction as F
import mzv_exact as M
import engine as E

# ---------- 2-var poly ops: dict (e1,e2)->F ----------
def padd(p,q,s=F(1)):
    out=dict(p)
    for e,c in q.items():
        out[e]=out.get(e,F(0))+s*c
        if out[e]==0: del out[e]
    return out
def pscale(p,s):
    return {} if s==0 else {e:c*s for e,c in p.items()}
def pmul(p,q):
    out={}
    for (a,b),c1 in p.items():
        for (d,e2),c2 in q.items():
            k=(a+d,b+e2); out[k]=out.get(k,F(0))+c1*c2
    return out
def linpow2(al,be,m):  # (al*x+be*y)^m
    from math import comb
    return {(k,m-k): F(comb(m,k))*al**k*be**(m-k) for k in range(m+1) if al**k*be**(m-k)}

# ---------- zmon series: dict zmon(sorted tuple)->poly ----------
def zs_add(a,b,s=F(1)):
    out=dict(a)
    for z,p in b.items():
        out[z]=padd(out.get(z,{}),p,s)
        if not out[z]: del out[z]
    return out
def zs_mul(a,b):
    out={}
    for z1,p1 in a.items():
        for z2,p2 in b.items():
            z=tuple(sorted(z1+z2))
            out[z]=padd(out.get(z,{}),pmul(p1,p2))
    return {z:p for z,p in out.items() if p}

def gamma_product_series(terms, Wmax):
    """terms: list of ((al,be), e) for prod Gamma(1+(al*x+be*y)*u)^e. Returns series E[W]."""
    q=[{}]*(Wmax+1)
    for m in range(2,Wmax+1):
        acc={}
        for (al,be),e in terms:
            acc=padd(acc, pscale(linpow2(al,be,m), e))
        q[m]={(m,): pscale(acc, F((-1)**m, m))}
    E=[{} for _ in range(Wmax+1)]
    E[0]={(): {(0,0):F(1)}}
    for W in range(1,Wmax+1):
        acc={}
        for i in range(2,W+1):
            acc=zs_add(acc, zs_mul(q[i], E[W-i]), F(i,W))
        E[W]=acc
    return E

# Plane definitions: B and G0 terms
planes={
 'Pi1': {'B': [(((1,0),1)),(((1,0),1)),(((0,1),1)),(((0,1),1)),(((2,0),-1)),(((0,2),-1))],
         'G0':[(((1,1),1)),(((1,1),1)),(((2,2),-1))]},
 'Pi2': {'B': [(((1,0),1)),(((0,1),1)),(((1,0),1)),(((0,1),1)),(((1,1),-1)),(((1,1),-1))],
         'G0':[(((1,1),1)),(((1,1),1)),(((2,2),-1))]},
}
WMAX=12
TABLES={}
def table(W):
    if W not in TABLES:
        TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return TABLES[W]

result={}
t0=time.time()
for name in ['Pi1','Pi2']:
    Bs=gamma_product_series(planes[name]['B'],WMAX)
    Gs=gamma_product_series(planes[name]['G0'],WMAX)
    H=[zs_add(Bs[W],Gs[W],F(1)) for W in range(WMAX+1)]
    H=[{z:pscale(p,F(1,2)) for z,p in h.items()} for h in H]
    # reduce to DS basis
    coords={}
    for W in range(2,WMAX+1):
        T=table(W); d=T.nullity
        acc=[{} for _ in range(d)]
        for zmon,poly in H[W].items():
            mon=tuple((k,) for k in zmon)
            wc=M.wordcombo(mon)
            for w,c in wc.items():
                rw=T.reduce_word(w)
                for b in range(d):
                    if rw[b]: acc[b]=padd(acc[b], pscale(poly,c*rw[b]))
        coords[W]=acc
    result[name]=coords
    print(f"{name}: series+reduction done, t={time.time()-t0:.1f}s")
pickle.dump(result, open('partB_planecoords.pkl','wb'))
print("saved. total", time.time()-t0)
