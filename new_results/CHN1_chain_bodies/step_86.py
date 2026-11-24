import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
from fractions import Fraction as F
import pickle, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(pp):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in pp.items())
P8=to_sp(forms[(8,3)])
gradP8=[sp.diff(P8,v).subs({x:1,y:1,z:1,w:1}) for v in (x,y,z,w)]
print('grad P8(1,1,1,1) =', gradP8, '=', sp.Rational(54,5), '*(1,-1,-1,1)?', [sp.Rational(g)/(sp.Rational(54,5)) for g in gradP8])
Q5=sp.cancel(sp.Rational(224)*to_sp({tuple(xx-1 for xx in e):c for e,c in forms[(10,6)].items()})/(x+y+z+w))
# print Q5 grouped by total-degree in a readable way
terms=[]
for mon,c in sorted(sp.Poly(Q5).terms(), key=lambda mc:(-sum(mc[0]), mc[0])):
    m='*'.join(v**e for v,e in zip((x,y,z,w),mon) if e>0)
    terms.append('%s%s'%(c if c<0 else '+'+str(c), m if m else ''))
print('Q5 =')
print(' '.join(terms))
