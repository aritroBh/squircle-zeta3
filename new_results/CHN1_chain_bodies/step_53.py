import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, symengine as S
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
lbl={10:{5:'z2*z35',6:'z37'},11:{7:'z3*z35',8:'z335'}}
print('=== coefficient VALUES (form / 1) at key points ===')
for pt in [(1,1,1,2),(3,1,1,1),(1,2,2,1)]:
    print(pt)
    for W,dd in lbl.items():
        for b,nm in dd.items():
            v=S.peval(forms[(W,b)], pt)
            print('  W%d %-7s = %s'%(W,nm,v))
# point on C3=0 with y=z=1, x=2, solve for w>0
C3 = to_sp({tuple(xx-1 for xx in e):c for e,c in forms[(8,3)].items()})/(x+y+z+w)*10
f=sp.Poly(C3.subs({x:2,y:1,z:1}),w)
roots=sp.nroots(f)
print('roots C3(2,1,1,w):', roots)
for rt in roots:
    if sp.im(rt)==0 and sp.re(rt)>0:
        pt=(2,1,1,float(sp.re(rt)))
        print('point', pt, ' on plane a1+a4=a2+a3?', abs(2+pt[3]-2)<1e-12)
        # evaluate z37 form and z35 form numerically
        for key,nm in [((8,3),'z35'),((10,6),'z37'),((10,5),'z2z35'),((11,8),'z335')]:
            v=sum(float(c)*pt[0]**e[0]*pt[1]**e[1]*pt[2]**e[2]*pt[3]**e[3] for e,c in forms[key].items())
            print('  %s form = %.3e'%(nm,v))
print('t=%.1fs'%(time.time()-t0))
