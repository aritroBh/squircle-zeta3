
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
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j) for j in (1,2,3)})
def tail2(poly):
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j)/2**(j+1) for j in (1,2,3)})
A2 = 2*z*(1-z); A3 = -4*c*z*(1-z); omc = 2*z
J1 = sp.expand(omc*avg(u1))
J2 = sp.expand(omc*avg(u2) - A2*tail1(u1)/a)
J3 = sp.expand(omc*avg(u3) - A2*(tail1(u2) - tail1(u1))/a - A3*tail2(u1)/a**2)

# isolate the (Xi,k)-independent part of J3
J3c = sp.expand(J3.subs({Xi:0, k2:0, k3:0}))
print("J3 const part:", sp.factor(J3c))
print("expect: -4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 =", sp.factor(-4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3))
print("match:", sp.simplify(J3c - (-4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3))==0)

# Tail functionals explicitly
print("\nT1(u1) =", tail1(u1), " T1(u2) =", tail1(u2), " T2(u1) =", tail2(u1))
print("-A2/a (T1u2-T1u1) =", sp.simplify(-A2*(tail1(u2)-tail1(u1))/a))
print("-A3/a^2 T2(u1) =", sp.simplify(-A3*tail2(u1)/a**2))

# Full pretty forms
print("\nJ2 =", sp.factor(J2))
J3_avg_part = sp.expand((J3 - J3c))
print("\nJ3 (averaged part) =", sp.factor(J3_avg_part.subs({k2: sp.Symbol('k2'), k3: sp.Symbol('k3')})))
