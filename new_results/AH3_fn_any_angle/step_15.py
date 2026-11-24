
import sympy as sp
z, c, a, A, B, Xi, k2, k3 = sp.symbols('z c a A B Xi k2 k3')
r = sp.Symbol('r')
d2 = (a**2 - A)/2
d3 = a**3/3 - a*A/2 + B/6
e1 = -A/a
e2 = B/(2*a) - A**2/(2*a**2)
u1 = r/a
u2 = sp.expand((e1*u1 - d2*u1**2)/a)
u3 = sp.expand((e1*u2 + e2*u1**2 - 2*d2*u1*u2 - d3*u1**3)/a)
mu = {1: Xi, 2: Xi**2 + k2, 3: Xi**3 + 3*k2*Xi + k3}
def avg(poly):
    return sp.expand(poly).subs({r**j: mu[j] for j in (1,2,3)}).subs(r, mu[1])
def tail1(poly):
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j) for j in (1,2,3)}).subs(r, 1)
def tail2(poly):
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j)/2**(j+1) for j in (1,2,3)}).subs(r, sp.Rational(1,4))
A2 = 2*z*(1-z); A3 = -4*c*z*(1-z); omc = 2*z
J1 = sp.expand(omc*avg(u1))
J2 = sp.expand(omc*avg(u2) - A2*tail1(u1)/a)
J3 = sp.expand(omc*avg(u3) - A2*(tail1(u2) - tail1(u1))/a - A3*tail2(u1)/a**2)

# Check 1: f(x)=x (a=1,A=B=0) at general angle vs the EXACT Ahmed-series expansion
# exact: 2-L = 2z(1-exp(-eps*Xi + a2 eps^2 + a3 eps^3)) + eps*Lam1/(n(1-eps)) + eps*Lam2/((2-eps)n^2) + ...
# with eps=1/(n-1), 1/n = eps/(1+eps). Through eps^3, in terms of Xi=log n + k1, mu_j:
k1 = sp.Symbol('k1')
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
exact_J1 = 2*z*mu1
exact_J2 = -z*mu2 + 2*z*(1-z)
exact_J3 = z*mu3/3 - c*z*(1-z)
sub = {a:1, A:0, B:0}
print("f=x: J1 match:", sp.simplify(J1.subs(sub) - exact_J1) == 0)
print("f=x: J2 match:", sp.simplify(J2.subs(sub) - exact_J2) == 0)
print("f=x: J3 match:", sp.simplify(J3.subs(sub) - exact_J3) == 0)

# Check 2 (right angle): J2, J3 for general f
subr = {z: sp.Rational(1,2), c: 0}
print("\nRight angle J2 =", sp.simplify(J2.subs(subr)))
print("Right angle J3 =", sp.simplify(J3.subs(subr)))

# store symbolic expressions for later numeric use
import pickle
pickle.dump({'J1':J1,'J2':J2,'J3':J3,'u1':u1,'u2':u2,'u3':u3}, open('Js.pkl','wb'))
print("\nsaved")
