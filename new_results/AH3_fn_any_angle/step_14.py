
import sympy as sp

# Symbols: angle theta -> z = sin^2(th/2), c = cos th; f data a=f'(1), A=f''(1), B=f'''(1)
z, c, a, A, B, Xi, k2, k3 = sp.symbols('z c a A B Xi k2 k3', positive=True)
r = sp.Symbol('r')  # \hat rho = log(n a / s)

# Level-set recursion: solve lambda(t) - eps log f'(1-t) = eps*rho, t = eps u1 + eps^2 u2 + eps^3 u3
d2 = (a**2 - A)/2
d3 = a**3/3 - a*A/2 + B/6
e1 = -A/a
e2 = B/(2*a) - A**2/(2*a**2)

u1 = r/a
u2 = (e1*u1 - d2*u1**2)/a
u3 = (e1*u2 + e2*u1**2 - 2*d2*u1*u2 - d3*u1**3)/a
u2 = sp.expand(u2); u3 = sp.expand(u3)
print("u1 =", u1)
print("u2 =", u2)
print("u3 =", u3)

# moments (cumulant form): mu1 = Xi, mu2 = Xi^2 + k2, mu3 = Xi^3 + 3 k2 Xi + k3
mu = {1: Xi, 2: Xi**2 + k2, 3: Xi**3 + 3*k2*Xi + k3}
def avg(poly):   # replace r^j by mu_j
    return sp.expand(poly).subs({r**j: mu[j] for j in (1,2,3)}).subs(r, mu[1])
def tail1(poly): # r^j -> (-1)^j j!   (1/n tail)
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j) for j in (1,2,3)}).subs(r, 1)
def tail2(poly): # r^j -> (-1)^j j!/2^{j+1}  (1/n^2 tail)
    return sp.expand(poly).subs({r**j: (-1)**j*sp.factorial(j)/2**(j+1) for j in (1,2,3)}).subs(r, sp.Rational(1,4))

A2 = 2*z*(1-z)          # = sin^2(theta)/2
A3 = -4*c*z*(1-z)       # = -cos(theta) sin^2(theta)
omc = 2*z               # 1 - cos(theta)

J1 = sp.expand(omc*avg(u1))
J2 = sp.expand(omc*avg(u2) - A2*tail1(u1)/1)   # note: -A2/a * tail1(u1); tail1(u1) includes 1/a
J3 = sp.expand(omc*avg(u3) - A2*(tail1(u2) - tail1(u1)) - A3*tail2(u1))
# careful: tail functionals already include the 1/a factors inside u_k; the outer factors:
# J2 = omc*avg(u2) - (A2)*tail1(u1)          [since eps*tail/n with 1/(na) absorbed]
print("\nJ1 =", J1)
print("J2 =", sp.expand(J2))
print("J3 =", sp.expand(J3))
