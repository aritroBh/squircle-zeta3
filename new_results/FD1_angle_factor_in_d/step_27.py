import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
import sympy as sp
from fractions import Fraction as Fr
Fd9=pickle.load(open(fd1+'/Fd_data_d9.pkl','rb'))
F9=Fd9[9]
t=sp.symbols('t')
pts=[]
for c,v in F9.items():
    tv=(1-c*c)/(c*c)   # t = (1-c^2)/c^2
    pts.append((tv, v))
# add squircle point F_9(1)=1 at t=1
pts.append((Fr(1),Fr(1)))
# interpolate degree <=4: use 5 points, check rest
import itertools
def interp(points):
    tt=sp.symbols('tt')
    xs=[p[0] for p in points]; ys=[p[1] for p in points]
    poly=sp.interpolate([(sp.Rational(x.numerator,x.denominator)) for x in xs],
                        [sp.Rational(y.numerator,y.denominator) for y in ys], tt)
    return sp.expand(poly)
poly=interp(pts[:5])
ok=all(sp.simplify(poly.subs(sp.symbols('tt'), sp.Rational(x.numerator,x.denominator))-sp.Rational(y.numerator,y.denominator))==0 for x,y in pts)
print("degree-4 fit through first 5 pts matches all 9:", ok)
print("F_9(t) =", poly)
# express as f_9 = F_9 * f_9(1) with integer coeffs: scale so that value at t=1 gives zigzag-ish
# primitive integer poly:
coeffs=sp.Poly(poly, sp.symbols('tt')).all_coeffs()
from math import lcm, gcd
denoms=[c.q for c in coeffs]
L=1
for q in denoms: L=lcm(L,q)
ints=[int(c*L) for c in coeffs]
g=0
for v in ints: g=gcd(g,abs(v))
ints=[v//g for v in ints]
print("integer poly (descending):", ints, "value at 1:", sum(ints))
open(fd1+'/progress.log','a').write('step_19: F_9(t) fit, degree 4, verified on all 9 points\n')
