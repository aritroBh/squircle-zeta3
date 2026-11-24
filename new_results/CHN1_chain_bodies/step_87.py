import os
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
from fractions import Fraction as F
import pickle, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(pp):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in pp.items())
Q5=sp.cancel(sp.Rational(224)*to_sp({tuple(xx-1 for xx in e):c for e,c in forms[(10,6)].items()})/(x+y+z+w))
vs='xyzw'
terms=[]
for mon,c in sorted(sp.Poly(Q5).terms(), key=lambda mc:(-sum(mc[0]), mc[0])):
    m=''.join(vs[i]+('' if mon[i]==1 else '^%d'%mon[i]) for i in range(4) if mon[i]>0)
    terms.append('%+s*%s'%(c,m))
print('Q5 =')
print(' '.join(terms))
