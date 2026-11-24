import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
open('progress.log','a').write('step_38-39: BLIND CHECK PASS: c8 at (3,1,1,1) matches prediction to 4.23e-71 (<1e-60), imag -2.5e-76, c0=1\n')
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, symengine as S
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w,s,t=sp.symbols('x y z w s t')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
def vantests(Q):
    qs=to_sp(Q)
    res={}
    for name,sub in {'Pi1':{x:s,y:s,z:t,w:t},'Pi2':{x:s,y:t,z:s,w:t},
                     'diag':{x:1,y:1,z:1,w:1},'a1+a4=a2+a3':{w:y+z-x},'pal':{x:s,y:t,z:t,w:s},
                     'a1=a2=a3':{x:s,y:s,z:s,w:t}}.items():
        v=sp.expand(qs.subs(sub)) if not isinstance(sub[x],int) else qs.subs(sub)
        res[name]=v
    return res
# zeta(3,7) form, divided by a1a2a3a4
P=forms[(10,6)]
Q={tuple(xx-1 for xx in e):c for e,c in P.items()}
res=vantests(Q)
print('=== zeta(3,7) form / (a1a2a3a4): degree 6 ===')
for k,v in res.items():
    print(('VANISHES ' if v==0 else 'NONZERO  ')+k, '' if v==0 else sp.factor(v))
print('t=%.1fs'%(time.time()-t0))
