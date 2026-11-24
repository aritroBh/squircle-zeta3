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
def restrict(W, sub):
    T=E.table(W)
    return [sp.expand(to_sp(forms[(W,b)]).subs(sub)) for b in range(T.nullity)]
def logcoeffs(clist, Wmax):
    logL={}
    for W in range(2,Wmax+1):
        acc = list(clist[W])
        for i in range(2,W-1):
            if W-i<2: continue
            acc = [a - sp.Rational(i,W)*li*ci for a,li,ci in zip(acc, logL[i], clist[W-i])]
        logL[W]=[sp.expand(a) for a in acc]
    return logL
labels={2:['z2'],3:['z3'],4:['z4'],5:['z5','z2z3'],6:['z6','z3^2'],7:['z7','z2z5','z4z3'],
        8:['z8','z3z5','z2z3^2','z35'],9:['z9','z2z7','z4z5','z3^3','z2^3z3'],
        10:['z10','z3z7','z5^2','z2z3z5','z4z3^2','z2z35','z37'],
        11:['z11','z2z9','z2^2z7','z2^3z5','z2^4z3','z2z3^3','z3^2z5','z3z35','z335']}
for name,sub in {'diag':{x:s,y:s,z:s,w:s}, 'Pi1':{x:s,y:s,z:t,w:t}, 'Pi2':{x:s,y:t,z:s,w:t}}.items():
    clist={W:restrict(W,sub) for W in range(2,12)}
    logL=logcoeffs(clist,11)
    print('=== locus',name,'===')
    for W in range(2,12):
        nz=[(labels[W][b], sp.factor(logL[W][b])) for b in range(len(logL[W])) if logL[W][b]!=0]
        only = len(nz)==1 and nz[0][0]=='z%d'%W
        print(' W=%d %s %s'%(W, 'SINGLE' if only else 'PRODUCTS', nz))
print('t=%.1fs'%(time.time()-t0))
