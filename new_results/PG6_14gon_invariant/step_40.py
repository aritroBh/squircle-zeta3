import pickle, sys, time
from fractions import Fraction as Fr
from math import comb, factorial
sys.path.insert(0,'.')
import ftfield as F
t0=time.time()
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
constvec = pickle.load(open('constvec14.pkl','rb'))
TERMS = ['a3','a4','a1','a2','a5','a6','a0']
WSLOT = {'a3':0,'a4':1,'a1':2,'a2':3,'a5':4,'a6':5,'a0':None}
# letter polynomial (LP): dict (f0,f1,f2)->Fr
def lp_mul(P,Q):
    R={}
    for (a,b,c),x in P.items():
        for (d,e,f),y in Q.items():
            k=(a+d,b+e,c+f); R[k]=R.get(k,Fr(0))+x*y
    return {k:v for k,v in R.items() if v}
def lp_scale(P,s):
    return {k:v*s for k,v in P.items() if v*s}
def lp_pow_lin(lin, j):
    # lin = (c2,ca,cu); (c2 l2+ ca la+ cu lu)^j
    if j==0: return {(0,0,0):Fr(1)}
    c2,ca,cu = lin
    R={}
    for f0 in range(j+1):
        for f1 in range(j-f0+1):
            f2 = j-f0-f1
            coef = Fr(factorial(j), factorial(f0)*factorial(f1)*factorial(f2)) * c2**f0 * ca**f1 * cu**f2
            if coef: R[(f0,f1,f2)]=coef
    return R
# precompute letter powers for each term
LP = {term:[lp_pow_lin(constvec[term], j) for j in range(4)] for term in TERMS}
mons14={}
def acc(mk, cFT):
    mons14[mk] = F.fadd(mons14.get(mk,F.ZERO), cFT)
for key7, poly in E3.items():
    # cur: dict (LP_dict, a,b,wc4) -> Fr coefficient on the letter-poly side is folded in; we carry (lpdict,a,b,wc)->Fr? 
    # Better: cur = dict (a,b,wc4) -> LP (letter poly with rational coeffs). Start: {(0,0,(0,0,0,0)): {(0,0,0):1}}
    cur = { (0,0,(0,0,0,0)): {(0,0,0):Fr(1)} }
    for ti,term in enumerate(TERMS):
        e=key7[ti]
        if e==0: continue
        slot=WSLOT[term]
        nxt={}
        for (a,b,wc),lp in cur.items():
            for k in range(e+1):
                if slot is None:
                    if k!=0: continue
                    na,nb,nwc = a,b,wc
                else:
                    na,nb,nwc = a,b,list(wc)
                    if slot==0: na=a+k
                    elif slot==1: nb=b+k
                    else: nwc[slot-2]+=k
                    nwc=tuple(nwc)
                add = lp_scale( lp_mul(lp, LP[term][e-k]), Fr(comb(e,k)) )
                key=(na,nb,nwc)
                if key not in nxt: nxt[key]={}
                for fk,fv in add.items():
                    nxt[key][fk]=nxt[key].get(fk,Fr(0))+fv
        cur = {k:{fk:fv for fk,fv in v.items() if fv} for k,v in nxt.items()}
        cur = {k:v for k,v in cur.items() if v}
    # multiply by t-poly [c0,c1,c2,c3] (FT) and accumulate
    for (a,b,wc),lp in cur.items():
        for fk,fv in lp.items():
            for n0 in range(4):
                pn=poly[n0]
                if pn[:3]==(Fr(0),)*3: continue
                acc((n0,(a,b),wc,fk), F.fscale(pn, fv))
print("mons14 size:", len(mons14), " time %.1fs"%(time.time()-t0))
# weight dist: a+b+sum(wc)+ (f0+f1+f2)
wd={}
for (n,ab,wc,fk) in mons14:
    w=ab[0]+ab[1]+sum(wc)+sum(fk); wd[w]=wd.get(w,0)+1
print("weight dist:", wd)
pickle.dump(mons14, open('mons14.pkl','wb'))
open('progress.log','a').write('PG6 s_B3: mons14 rebuilt with letter polys, %d monomials, weights %s\n'%(len(mons14),wd))
