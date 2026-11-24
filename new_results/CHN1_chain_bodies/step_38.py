import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, symengine as S, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
def info(name, P):
    div=all(all(x_>=1 for x_ in e) for e in P)
    print('%s: #monomials=%d, divisible by a1a2a3a4: %s'%(name,len(P),div))
    if div:
        Q={tuple(xx-1 for xx in e):c for e,c in P.items()}
        qs=to_sp(Q)
        t=time.time(); f=sp.factor(qs)
        print('  reduced degree-%d form factors as:'%(sum(next(iter(Q))) if Q else 0), f, '(%.0fs)'%(time.time()-t))
        return Q,f
    else:
        qs=to_sp(P); t=time.time(); f=sp.factor(qs)
        print('  factors as:', f, '(%.0fs)'%(time.time()-t))
        return None,f
labels={10:{5:'z2*z35',6:'z37'},11:{7:'z3*z35',8:'z335'}}
out={}
for W,dd in labels.items():
    for b,nm in dd.items():
        Q,f=info('W%d %s'%(W,nm), forms[(W,b)])
        out[(W,b)]=(Q,f)
print('t=%.1fs'%(time.time()-t0))
