import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, engine as E, mzv_exact as M
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w,s,t=sp.symbols('x y z w s t')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
# value vectors: W -> dict monomial(tuple of letters) -> sympy coeff, restricted to locus
def value_vec(W, sub):
    T=E.table(W)
    out={}
    for b,mon in enumerate(T.monomials):
        cf=sp.expand(to_sp(forms[(W,b)]).subs(sub))
        if cf!=0: out[mon]=cf
    return out
def val_mul(v1,v2):
    out={}
    for m1,c1 in v1.items():
        for m2,c2 in v2.items():
            m=tuple(sorted(m1+m2))
            out[m]=out.get(m,0)+c1*c2
    return out
def val_scale(v,c): return {m:cc*c for m,cc in v.items()}
def val_add(v1,v2):
    out=dict(v1)
    for m,c in v2.items(): out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items() if c!=0}
def val_reduce(v, W):
    # reduce monomial->DS coords via wordcombo; returns list of sympy exprs
    T=E.table(W)
    acc=[sp.Integer(0)]*T.nullity
    for mon,cf in v.items():
        wc=M.wordcombo(mon)
        coords=T.reduce(wc)
        for b in range(T.nullity):
            if coords[b]: acc[b]=acc[b]+cf*sp.Rational(coords[b].numerator,coords[b].denominator)
    return [sp.expand(a) for a in acc]
labels={2:['z2'],3:['z3'],4:['z4'],5:['z5','z2z3'],6:['z6','z3^2'],7:['z7','z2z5','z4z3'],
        8:['z8','z3z5','z2z3^2','z35'],9:['z9','z2z7','z4z5','z3^3','z2^3z3'],
        10:['z10','z3z7','z5^2','z2z3z5','z4z3^2','z2z35','z37'],
        11:['z11','z2z9','z2^2z7','z2^3z5','z2^4z3','z2z3^3','z3^2z5','z3z35','z335']}
for name,sub in {'diag':{x:s,y:s,z:s,w:s}, 'Pi1':{x:s,y:s,z:t,w:t}, 'Pi2':{x:s,y:t,z:s,w:t}}.items():
    cv={W:value_vec(W,sub) for W in range(2,12)}
    logV={}
    print('=== locus',name,'===')
    for W in range(2,12):
        acc=dict(cv[W])
        for i in range(2,W-1):
            if W-i<2: continue
            acc=val_add(acc, val_scale(val_mul(logV[i],cv[W-i]), -sp.Rational(i,W)))
        logV[W]={m:c for m,c in acc.items() if c!=0}
        red=val_reduce(logV[W],W)
        nz=[(labels[W][b], sp.factor(red[b])) for b in range(len(red)) if red[b]!=0]
        only = len(nz)==1 and nz[0][0]=='z%d'%W
        print(' W=%d %s %s'%(W, 'SINGLE' if only else 'PRODUCTS', nz))
print('t=%.1fs'%(time.time()-t0))
