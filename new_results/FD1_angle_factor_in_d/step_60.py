import sympy as sp
from fractions import Fraction as Fr
x=sp.symbols('x')
f={2:[0,1],3:[1,3],4:[0,2,3],5:[2,15,15],6:[0,17,60,45],7:[17,231,525,315],
   8:[0,62,378,630,315],9:[62,1320,5040,6615,2835],10:[0,1382,12720,34965,37800,14175]}
A={2:1,3:2,4:5,5:16,6:61,7:272,8:1385,9:7936,10:50521}
P={}
for d in range(2,11):
    m=d//2
    poly=sum(Fr(a)*(1-x)**j*(1+x)**(m-j) for j,a in enumerate(f[d]))
    P[d]=sp.expand(poly*Fr(A[d],sum(f[d])))
    print(f"P_{d} =", P[d], " P(-1)=",P[d].subs(x,-1)," P(1)=",P[d].subs(x,1))
# check P_d(-1) = d!/2^floor((d-1)/2)
import math
for d in range(2,11):
    print(d, P[d].subs(x,-1), math.factorial(d)//2**((d-1)//2))
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
open(fd1+'/progress.log','a').write('step_57: P_d table through d=10; P_10(-1) matches d!/2^4 pattern\n')
