import sympy as sp
i=sp.I; sq2=sp.sqrt(2); a=1+sq2; b=sq2-1; v=sp.symbols('v')
# fact 1: 1/(1+a^2)+1/(1+b^2)=1
print("1/(1+a^2)+1/(1+b^2) =", sp.simplify(1/(1+a**2)+1/(1+b**2)))
# fact 2: -2[(a^2-1)/(a^2+1)^2+(b^2-1)/(b^2+1)^2] = -1
print("-2[(a2-1)/(a2+1)^2+(b2-1)/(b2+1)^2] =", sp.simplify(-2*((a**2-1)/(a**2+1)**2+(b**2-1)/(b**2+1)**2)))
# fact 3: t'(i)=0, t''(i)=-i/2
t=(v**2+2*v-1)/(4*v)
print("t'(i)=", sp.simplify(sp.diff(t,v).subs(v,i)), " t''(i)=", sp.simplify(sp.diff(t,v,2).subs(v,i)))
# fact 4: K_{d+1}/K_d = -i(2d-1)/(d+1)   (I1, used everywhere)
from sympy import factorial, factorial2
d=sp.symbols('d', integer=True, positive=True)
def K(D): return (-1)**D*factorial2(2*D-3)*i**(D+1)/(2*sq2*factorial(D))
# numeric check for many d
ok=True
for D in range(2,30):
    lhs = sp.simplify(K(D+1)/K(D)); rhs = -i*(2*D-1)/(D+1)
    if sp.simplify(lhs-rhs)!=0: ok=False; print("I1 fail",D)
print("I1 K_{d+1}/K_d=-i(2d-1)/(d+1), d=2..29:", ok)
# algebraic sum check: (2d-1-D)+D+(1-d)/2-1 = (3d-3)/2
Dsym=sp.symbols('D')
print("sum check:", sp.simplify((2*d-1-Dsym)+Dsym+(1-d)/2-1 - (3*d-3)/2)==0)
with open("<project>/research/08_theory/TOW1R_review/progress.log","a") as fh:
    fh.write("step_20: rho2 derivation algebraic pieces (1/(1+a2)+1/(1+b2)=1, t'(i)=0,t''(i)=-i/2, I1, sum) verified\n")
