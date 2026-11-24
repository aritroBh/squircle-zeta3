import sympy as sp
y,x=sp.symbols('y x')
N=14
# p(d,1) values d=2..9
p1={2:1,3:1,4:6,5:13,6:89,7:297,8:2396,9:10841}
import math
H=sum(sp.Rational(p1[d], math.factorial(d))*y**(d-2) for d in range(2,10))
sy=sp.series(sp.sin(y),y,0,N+2).removeO().expand()
cy=sp.series(sp.cos(y),y,0,N+2).removeO().expand()
secy=sp.series(1/cy,y,0,N+2).removeO().expand()
h=sp.expand(H*(cy**2)/(1+sy))
h=sp.series(h,y,0,N).removeO().expand()
print("h(y) =", h)
# check k=2: p(d,2) = d! [y^{d-2k}] (1+sin y) h^k sec^{k+1} y
p2={4:1,5:1,6:29,7:60,8:1134,9:3651}
pred2=sp.expand((1+sy)*h**2*secy**3)
for d in range(4,10):
    coef=sp.expand(pred2).coeff(y,d-4)*math.factorial(d)
    print(f"d={d}: p(d,2) actual={p2[d]} pred={sp.nsimplify(coef)}")
