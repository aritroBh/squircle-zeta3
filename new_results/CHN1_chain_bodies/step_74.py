import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
Q4={tuple(xx-1 for xx in e):c for e,c in forms[(8,3)].items()}
C3=sp.cancel(10*to_sp(Q4)/(x+y+z+w))
P37=forms[(10,6)]; Q37={tuple(xx-1 for xx in e):c for e,c in P37.items()}
Q5=sp.cancel(sp.Rational(224)*to_sp(Q37)/(x+y+z+w))
grad=[sp.diff(Q5,v).subs({x:1,y:1,z:1,w:1}) for v in (x,y,z,w)]
print('grad Q5 at (1,1,1,1):', grad)
# C3 and Q5 on z=1, eliminate w
t=time.time()
C3z=sp.Poly(C3.subs({z:1}),x,y,w)
Q5z=sp.Poly(Q5.subs({z:1}),x,y,w)
print('degrees in w:', sp.degree(C3z,w), sp.degree(Q5z,w))
R=sp.resultant(C3z,Q5z,w)
print('resultant computed, total degree ~', sp.Poly(R).total_degree() if R!=0 else 0,'(%.0fs)'%(time.time()-t))
f=sp.factor(R)
print('factor(R) =', f)
print('t=%.1fs'%(time.time()-t0))
