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
found=[]
for xv in [sp.Rational(1,2),sp.Rational(1,3),sp.Rational(3),sp.Rational(4)]:
    for yv in [sp.Rational(1,2),sp.Rational(2),sp.Rational(1,3)]:
        for zv in [sp.Rational(1,2),sp.Rational(2),sp.Rational(3)]:
            f=sp.Poly(C3.subs({x:xv,y:yv,z:zv}),w)
            for rt in sp.nroots(f):
                if abs(sp.im(rt))<1e-30 and sp.re(rt)>1e-9:
                    wv=float(sp.re(rt))
                    if abs(float(xv)+wv-float(yv)-float(zv))>1e-6:  # off the plane a1+a4=a2+a3
                        found.append((float(xv),float(yv),float(zv),wv))
print('positive C3=0 points off the plane:', found[:6])
if found:
    pt=found[0]
    print('using', pt, ' plane check a1+a4-(a2+a3) =', pt[0]+pt[3]-pt[1]-pt[2])
    for key,nm in [((8,3),'z35'),((10,6),'z37'),((10,5),'z2z35'),((11,8),'z335'),((11,7),'z3z35')]:
        v=sum(float(c)*pt[0]**e[0]*pt[1]**e[1]*pt[2]**e[2]*pt[3]**e[3] for e,c in forms[key].items())
        print('  %-6s form = %.4e'%(nm,v))
print('t=%.1fs'%(time.time()-t0))
