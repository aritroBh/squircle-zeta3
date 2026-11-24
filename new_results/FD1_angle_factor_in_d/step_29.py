import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
import sympy as sp
x=sp.symbols('x')
# f_d integer polys in t (ascending coeffs)
f={2:[0,1],3:[1,3],4:[0,2,3],5:[2,15,15],6:[0,17,60,45],7:[17,231,525,315],
   8:[0,62,378,630,315],9:[62,1320,5040,6615,2835]}
A={2:1,3:2,4:5,5:16,6:61,7:272,8:1385,9:7936}
t=(1-x)/(1+x)
P={}
for d in range(2,10):
    m=d//2
    poly=sum(Fr(a)*(1-x)**j*(1+x)**(m-j) for j,a in enumerate(f[d]))
    poly=sp.expand(poly)
    # P_d = A_d/f_d(1) * poly
    f1=sum(f[d])
    P[d]=sp.expand(poly*Fr(A[d],f1))
    print(f"P_{d} =", P[d], "  P(0)=",P[d].subs(x,0)," P(-1)=",P[d].subs(x,-1)," P(1)=",P[d].subs(x,1))
open(fd1+'/progress.log','a').write('step_30: P_d(ca) table completed through d=9\n')
