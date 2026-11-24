# step_46: decagon canonical pair (fixed L(chi5,3))
import mpmath as mp
mp.mp.dps = 100
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
def chi5(n):
    r = n % 5
    return 0 if r==0 else (1 if r in (1,4) else -1)
L5 = mp.nsum(lambda n: chi5(int(n))/mp.mpf(n)**3, [1, mp.inf])
print("L(chi5,3) =", L5)
z3 = mp.zeta(3)
phi = (1+mp.sqrt(5))/2
X1, Xm1 = L3(1/phi), L3(-1/phi)
alpha = (X1+Xm1)/(2*z3)
beta  = (X1-Xm1)/(2*L5)
import sympy as sp
print("alpha =", alpha, "->", sp.nsimplify(str(alpha)))
print("beta  =", beta, "->", sp.nsimplify(str(beta)))
# local zeta(3) in S: C5/(60 L_edge), C5 = 48*5*sin^4(pi/10)*cos(pi/10)*(3 sin^2(pi/10)-1)
s10 = mp.sin(mp.pi/10); c10 = mp.cos(mp.pi/10)
C5 = 48*5*s10**4*c10*(3*s10**2-1)
Ledge = mp.sqrt((25-11*mp.sqrt(5))/10)
c3local = C5/(60*Ledge)
print("C5 =", C5, " c3local in S =", c3local, "->", sp.nsimplify(str(c3local), [sp.sqrt(5)]))
# N3 = (c1+c2)*alpha*zeta3 + (c1-c2)*beta*L5 + (c3 - c3local)*zeta3
from fractions import Fraction as Fr
s5 = sp.sqrt(5)
c1 = sp.Rational(-75,32) + sp.Rational(81,64)*s5
c2 = sp.Rational(-35,32) + sp.Rational(53,64)*s5
c3 = sp.Rational(33,64) - sp.Rational(19,80)*s5
a_rat = sp.Rational(sp.nsimplify(str(alpha)))
b_rat = sp.Rational(sp.nsimplify(str(beta)))
c3loc_rat = sp.nsimplify(str(c3local), [s5])
a = (c1+c2)*a_rat + c3 - c3loc_rat
b = (c1-c2)*b_rat
print("a (zeta(3) coefficient of N3) =", sp.simplify(a))
print("b (L(chi5,3) coefficient of N3) =", sp.simplify(b))
# numeric check: a*zeta3 + b*L5 vs direct combination
aval = float(sp.N(a, 30)); bval = float(sp.N(b, 30))
lhs = aval*z3 + bval*L5
rhs = (c1+c2)*a_rat*z3 + (c1-c2)*b_rat*L5 + (c3-c3loc_rat)*z3
rhs = mp.mpf(str(sp.N(rhs, 40)))
print("lhs-rhs:", lhs-rhs)
