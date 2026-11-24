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
cases=[('y=x+1', {y:x+1,z:1}, x), ('y=1-x', {y:1-x,z:1}, x), ('x=1 (Pi2)', {x:1,z:1}, y), ('y=x (Pi1)', {y:x,z:1}, x),
       ('x=-1', {x:-1,z:1}, y), ('x=-2', {x:-2,z:1}, y), ('y=-x', {y:-x,z:1}, x), ('y=-x/2', {y:-x/2,z:1}, x),
       ('x=-1/2', {x:sp.Rational(-1,2),z:1}, y), ('y=-2x', {y:-2*x,z:1}, x), ('y=-1', {y:-1,z:1}, x), ('y=-2x-1', {y:-2*x-1,z:1}, x),
       ('x=0', {x:0,z:1}, y), ('y=0', {y:0,z:1}, x)]
for name,sub,v1 in cases:
    c=sp.Poly(C3.subs(sub), v1, w); q=sp.Poly(Q5.subs(sub), v1, w)
    g=sp.gcd(c,q,v1,w)
    if sp.degree(g,w)>=1:
        sol=sp.solve(sp.Eq(g.as_expr(),0), w)
        print('%-12s component(s): w = %s'%(name, sol))
    else:
        print('%-12s (no common w-component)'%name)
print('t=%.1fs'%(time.time()-t0))
