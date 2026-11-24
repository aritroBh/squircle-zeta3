
import pickle, time
from fractions import Fraction as F
import mzv_exact as M

def padd(p,q,s=F(1)):
    out=dict(p)
    for e,c in q.items():
        out[e]=out.get(e,F(0))+s*c
        if out[e]==0: del out[e]
    return out
def pscale(p,s): return {} if s==0 else {e:c*s for e,c in p.items()}
def pmul(p,q):
    out={}
    for (a,b),c1 in p.items():
        for (d,e2),c2 in q.items():
            k=(a+d,b+e2); out[k]=out.get(k,F(0))+c1*c2
    return out
def linpow2(al,be,m):
    from math import comb
    return {(k,m-k): F(comb(m,k))*al**k*be**(m-k) for k in range(m+1) if al**k*be**(m-k)}
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

# CORRECTED planes
planes={
 'Pi1': {'B': [((1,0),1),((1,0),1),((0,1),1),((0,1),1),((2,0),-1),((0,2),-1)],
         'G0':[((1,1),1),((1,1),1),((2,2),-1)]},
 'Pi2': {'B': [((1,0),1),((0,1),1),((1,0),1),((0,1),1),((1,1),-1),((1,1),-1)],
         'G0':[((2,0),1),((0,2),1),((2,2),-1)]},
}
WMAX=12
TABLES={}
def table(W):
    if W not in TABLES:
        TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return TABLES[W]
result={}
for name in ['Pi1','Pi2']:
    Bs=gamma_product_series(planes[name]['B'],WMAX)
    Gs=gamma_product_series(planes[name]['G0'],WMAX)
    H=[{z:pscale(p,F(1,2)) for z,p in zs_add(Bs[W],Gs[W]).items()} for W in range(WMAX+1)]
    coords={}
    for W in range(2,WMAX+1):
        T=table(W); d=T.nullity
        acc=[{} for _ in range(d)]
        for zmon,poly in H[W].items():
            wc=M.wordcombo(tuple((k,) for k in zmon))
            for w,c in wc.items():
                rw=T.reduce_word(w)
                for b in range(d):
                    if rw[b]: acc[b]=padd(acc[b], pscale(poly,c*rw[b]))
        coords[W]=acc
    result[name]=coords
    print(name,"done")
pickle.dump(result, open('partB_planecoords.pkl','wb'))

# compare W=2..11 with forms.pkl
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
forms = pickle.load(open(f'{chn1}/forms.pkl','rb'))
def parse_form(v): return {tuple(int(x) for x in ks.split(',')):F(c) for ks,c in v.items()}
def restrict4(p4, plane):
    out={}
    for (e1,e2,e3,e4),c in p4.items():
        k=(e1+e2, e3+e4) if plane=='Pi1' else (e1+e3, e2+e4)
        out[k]=out.get(k,F(0))+c
    return {k:c for k,c in out.items() if c}
allok=True
for W in range(2,12):
    for plane in ['Pi1','Pi2']:
        d=len(result[plane][W]); bad=[]
        for b in range(d):
            fr=restrict4(parse_form(forms[(W,b)]), plane)
            diff=dict(fr)
            for e,c in result[plane][W][b].items():
                diff[e]=diff.get(e,F(0))-c
                if diff[e]==0: del diff[e]
            if diff: bad.append(b)
        if bad: allok=False
        print(f"W={W:>2} {plane}: {'OK' if not bad else 'MISMATCH '+str(bad)}")
print("ALL EXACT MATCH W=2..11:", allok)
