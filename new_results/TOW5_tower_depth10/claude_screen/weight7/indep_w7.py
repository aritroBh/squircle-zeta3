import pickle
from fractions import Fraction as F
W='<project>/research/08_theory/WTS_weight_k_structure/theory'
K={k:pickle.load(open(f'{W}/kern_{k}.pkl','rb')) for k in range(5,15) if True} if False else {}
import os
for k in range(1,19):
    p=f'{W}/kern_{k}.pkl'
    if os.path.exists(p): K[k]=pickle.load(open(p,'rb'))
print('levels',sorted(K))
Z=(F(0),)*4
def dif(seq,n):
    for _ in range(n): seq=[seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def test(d,kmin,kmax):
    keys=[(i,d-i) for i in range(d,(d-1)//2,-1)]
    out={}
    for key in keys:
        poles=set()
        for k in range(kmin,kmax+1): poles|=set(K[k][key]['real'])
        # pole presence uniform?
        pres={p:[p in K[k][key]['real'] for k in range(kmin,kmax+1)] for p in poles}
        nonuni=[p for p,v in pres.items() if len(set(v))>1]
        nd=0;maxdeg=-1;n=0
        for p in poles:
            for c in range(4):
                seq=[K[k][key]['real'].get(p,Z)[c] for k in range(kmin,kmax+1)]
                n+=1
                assert all(x==0 for x in dif(seq,d)),(key,p,c)
                for deg in range(d,-1,-1):
                    if any(x!=0 for x in dif(seq,deg)): maxdeg=max(maxdeg,deg);break
        ip=all(not K[k][key]['ipos'] and not K[k][key]['ineg'] for k in range(kmin,kmax+1))
        cst=[K[k][key]['const'] for k in range(kmin,kmax+1)]
        out[key]=(n,len(poles),len(nonuni),maxdeg,ip,set(map(str,cst)))
        print(d,key,'comps',n,'poles',len(poles),'nonuniform-presence',len(nonuni),'max nonzero diff order',maxdeg,'(=> degree',maxdeg-0 if maxdeg<0 else maxdeg,') iPoles empty',ip,'const',set(map(str,cst)))
print('--- d=6, k=7..13: d-th (6th) diffs zero; first nonzero diff order => degree+1')
test(6,7,13)
print('--- d=6 including k=6 : should FAIL')
try: test(6,6,13)
except AssertionError as e: print('k=6 included: 6th diff nonzero at',e)
# k=6 mismatches against quintic fit from 8..13 and also count against 7..12 fit
def lag(ys,xs,x):
    o=F(0)
    for a,xa in enumerate(xs):
        L=F(1)
        for b,xb in enumerate(xs):
            if b!=a: L*=F(x-xb,xa-xb)
        o+=ys[a]*L
    return o
for fitk in ([8,9,10,11,12,13],[7,8,9,10,11,12],[9,10,11,12,13,14]):
    if not all(k in K for k in fitk): continue
    for tgt in (6,7):
        if tgt in fitk: continue
        bad=tot=0
        for key in [(6,0),(5,1),(4,2),(3,3)]:
            poles=set()
            for k in fitk+[tgt]: poles|=set(K[k][key]['real'])
            for p in poles:
                for c in range(4):
                    tot+=1
                    if lag([K[k][key]['real'].get(p,Z)[c] for k in fitk],fitk,tgt)!=K[tgt][key]['real'].get(p,Z)[c]: bad+=1
        print('fit',fitk,'target',tgt,'bad',bad,'/',tot)
print('--- controls: d=5 (quartic) k=6..14 ; d=4 k=5..14 ; d=5 k=5 should fail')
test(5,6,14); test(4,5,14)
# does the k=7 d=6 kernel even exist for key (7,0)? degree-7 keys exist at k=7
print('keys at k=7 degree 7:',[k for k in K[7] if sum(k)==7])
print('degree of d=6 sequences: exact degree 5?')
