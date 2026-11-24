import pickle, sys, time
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0,'.')
import ftfield as F
t0=time.time()
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
constvec = pickle.load(open('constvec14.pkl','rb'))
TERMS = ['a3','a4','a1','a2','a5','a6','a0']
# wlog slot: 0=lt(a),1=l1t(b), 2..5 = wc[0..3]=(w1,w2,w5,w6); None=const
WSLOT = {'a3':0,'a4':1,'a1':2,'a2':3,'a5':4,'a6':5,'a0':None}
ZEROV=(0,0,0)
def vadd(u,v): return (u[0]+v[0],u[1]+v[1],u[2]+v[2])
def vscale(v,k): return (v[0]*k,v[1]*k,v[2]*k)
def faddto(poly,n,c): poly[n]=F.fadd(poly.get(n,F.ZERO),c)

mons = {}  # (letters, a, b, wc4) -> dict n->FT
for key7, poly in E3.items():
    # cur: (letters, wtuple6) -> FT coeff (scalar, before t-poly)
    cur = { (ZEROV,(0,0,0,0,0,0)): F.ONE }
    for ti,term in enumerate(TERMS):
        e = key7[ti]
        if e==0: continue
        lv = constvec[term]; slot = WSLOT[term]
        nxt = {}
        for (L,W),c0 in cur.items():
            if slot is None:
                L2 = vadd(L, vscale(lv,e))
                k2=(L2,W); nxt[k2]=F.fadd(nxt.get(k2,F.ZERO),c0)
            else:
                for k in range(e+1):
                    L2 = vadd(L, vscale(lv,e-k))
                    W2 = list(W); W2[slot]+=k; W2=tuple(W2)
                    cc = F.fscale(c0, Fr(comb(e,k)))
                    k2=(L2,W2); nxt[k2]=F.fadd(nxt.get(k2,F.ZERO),cc)
        cur = nxt
    # now multiply by t-poly [c0,c1,c2,c3] and store
    for (L,W),c in cur.items():
        a=W[0]; b=W[1]; wc=(W[2],W[3],W[4],W[5])
        mk=(L,a,b,wc)
        if mk not in mons: mons[mk]={}
        for n0 in range(4):
            pn = poly[n0]
            if pn[:3]==(Fr(0),)*3: continue
            faddto(mons[mk], n0, F.fmul(c,pn))
# flatten to mons14: (n,(a,b),wc,letters)->FT, keeping only nonzero
mons14={}
for (L,a,b,wc),poly in mons.items():
    for n,c in poly.items():
        if c[:3]==(Fr(0),)*3: continue
        mons14[(n,(a,b),wc,L)]=c
print("mons14 size:", len(mons14), " time %.1fs"%(time.time()-t0))
# weight distribution: a+b+sum(wc)+sum(|letters|)  (letter degree = sum of exponents, but letters are positive? they're the L vec which can be negative exponents -> use sum)
wd = {}
for (n,ab,wc,L) in mons14:
    w = ab[0]+ab[1]+sum(wc)+sum(L)
    wd[w]=wd.get(w,0)+1
print("weight (a+b+wc+lettersum) distribution:", wd)
pickle.dump(mons14, open('mons14.pkl','wb'))
open('progress.log','a').write('PG6 s_B2: mons14 built, %d reduced monomials, weights %s\n'%(len(mons14),wd))
