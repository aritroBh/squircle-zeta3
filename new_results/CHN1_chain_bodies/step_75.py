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
for name, ys in [('y=x+1', x+1), ('y=1-x', 1-x), ('x=1(Pi2)', y), ('y=x(Pi1)', x)]:
    c=sp.Poly(C3.subs({y:ys,z:1}), x, w) if name!='x=1(Pi2)' else sp.Poly(C3.subs({x:1,z:1}), y, w)
    q=sp.Poly(Q5.subs({y:ys,z:1}), x, w) if name!='x=1(Pi2)' else sp.Poly(Q5.subs({x:1,z:1}), y, w)
    v1 = y if name=='x=1(Pi2)' else x
    g=sp.gcd(c,q, v1, w)
    print('===', name, ': gcd degree in w:', sp.degree(g, w))
    if sp.degree(g,w)>=1:
        # solve gcd=0 for w as function of v1
        cw=sp.Poly(g, w)
        sol=sp.solve(sp.Eq(g,0), w)
        print('   component(s): w =', sol)
print('t=%.1fs'%(time.time()-t0))
