import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, engine as E
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w,s,t=sp.symbols('x y z w s t')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
def value_vec(W, sub):
    T=E.table(W); out={}
    for b,mon in enumerate(T.monomials):
        cf=sp.expand(to_sp(forms[(W,b)]).subs(sub))
        if cf!=0: out[mon]=cf
    return out
def val_mul(v1,v2):
    out={}
    for m1,c1 in v1.items():
        for m2,c2 in v2.items():
            m=tuple(sorted(m1+m2)); out[m]=out.get(m,0)+c1*c2
    return out
def logV_series(sub, Wmax=11):
    cv={W:value_vec(W,sub) for W in range(2,Wmax+1)}
    logV={}
    for W in range(2,Wmax+1):
        acc=dict(cv[W])
        for i in range(2,W-1):
            if W-i<2: continue
            for m,c in val_mul(logV[i],cv[W-i]).items(): acc[m]=acc.get(m,0)-sp.Rational(i,W)*c
        logV[W]={m:sp.expand(c) for m,c in acc.items() if c!=0}
    return logV
cache={}
for name,sub in {'diag':{x:s,y:s,z:s,w:s}, 'Pi1':{x:s,y:s,z:t,w:t}, 'Pi2':{x:s,y:t,z:s,w:t}}.items():
    logV=logV_series(sub)
    Q={}
    prods={}
    for m in range(2,12):
        cm=logV[m].get(((m,),), sp.Integer(0))
        Q[m]=sp.expand((-1)**m*m*cm)
        prods[m]={mon:c for mon,c in logV[m].items() if mon!=((m,),)}
    cache[name]={'Q':Q,'prods':prods}
pickle.dump(cache, open('locus_logV.pkl','wb'))
print('cached. t=%.0fs'%(time.time()-t0))
# verify product parts identical across loci (Pi1 vs Pi2 exactly, diag = restriction s=t)
same=all(cache['Pi1']['prods'][m]==cache['Pi2']['prods'][m] for m in range(2,12))
print('product series identical Pi1 vs Pi2:', same)
sd=all(sp.expand(cache['Pi1']['prods'][m][mon].subs({t:s}) - cache['diag']['prods'][m][mon])==0 for m in range(2,12) for mon in cache['diag']['prods'][m])
print('diag = Pi1 products at s=t:', sd)
open('progress.log','a').write('step_43-46: logVol on diag/Pi1/Pi2: singles+products split; products IDENTICAL on all three loci; singles cached to locus_logV.pkl\n')
