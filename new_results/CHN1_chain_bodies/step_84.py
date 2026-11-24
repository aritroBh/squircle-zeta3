import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(pp):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in pp.items())
Q4={tuple(xx-1 for xx in e):c for e,c in forms[(8,3)].items()}
C3=sp.cancel(10*to_sp(Q4)/(x+y+z+w))
Sa=x+y+z+w
print('z3z35 identity:', sp.expand(to_sp(forms[(11,7)]) - x*y*z*w*(y+w)*(x+z)*Sa**2*C3/10)==0)
mid = w**2*x + w*x**2 + 2*w*x*y + 2*w*x*z + 2*w*y*z + 2*x*y*z + y**2*z + y*z**2
print('z335 identity:', sp.expand(to_sp(forms[(11,8)]) + x*y*z*w*Sa*mid*C3/10)==0)
Q5=sp.cancel(sp.Rational(224)*to_sp({tuple(xx-1 for xx in e):c for e,c in forms[(10,6)].items()})/Sa)
print('z37 identity:', sp.expand(to_sp(forms[(10,6)]) - x*y*z*w*Sa*Q5/224)==0)
print('Q5 irreducible over Q:', sp.factor(Q5)==Q5)
print('Q5(1,1,1,1)=', Q5.subs({x:1,y:1,z:1,w:1}))
print('t=%.1fs'%(time.time()-t0))
