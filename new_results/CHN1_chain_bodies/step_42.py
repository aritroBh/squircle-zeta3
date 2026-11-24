import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
from fractions import Fraction as F
import pickle, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w,s,t=sp.symbols('x y z w s t')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
P=forms[(10,6)]
Q={tuple(xx-1 for xx in e):c for e,c in P.items()}
qs=to_sp(Q)
print('deg', sp.Poly(qs).total_degree() if qs else 0, '#terms', len(sp.Poly(qs).terms()) if qs else 0)
subd={'Pi1':{x:s,y:s,z:t,w:t},'Pi2':{x:s,y:t,z:s,w:t},'diag':{x:1,y:1,z:1,w:1},
      'a1+a4=a2+a3':{w:y+z-x},'pal':{x:s,y:t,z:t,w:s},'a1=a2=a3':{x:s,y:s,z:s,w:t}}
for name,sub in subd.items():
    v=sp.expand(qs.subs(sub))
    print(('VANISHES ' if v==0 else 'NONZERO  ')+name, '' if v==0 else sp.factor(v))
