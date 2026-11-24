
import sympy as sp

eps, x = sp.symbols('eps, x')
zs = sp.symbols('zs')
N = 1/eps**2
c3 = sp.Rational(-1,2)

# exact saddle
us = (-zs + sp.sqrt(zs**2 + 4*zs/eps**2))/2

# phi and its derivatives at us
u = sp.symbols('u')
L = sp.log(u+zs) - sp.log(u)
phi = -u - (1/eps**2)*L
dervs = [sp.diff(phi, u, j) for j in range(1,6)]
at = lambda expr, ord_: sp.series(sp.series(expr.subs(u, us), eps, 0, ord_).removeO().expand(), eps, 0, ord_)
phi_s = at(phi, 3)
phi2  = at(dervs[1], 4)
phi3  = at(dervs[2], 4)
phi4  = at(dervs[3], 4)
phi5  = at(dervs[4], 4)
print("phi_s =", phi_s)
print("phi2  =", phi2)
print("phi3  =", phi3)
print("phi4  =", phi4)

P = (u+zs)**sp.Rational(3,2)/u
Q = 1 + sp.Rational(3,2)*c3*(u+zs)*eps**2
F  = P*Q
F_s = at(F, 3)
F1  = at(sp.diff(F,u,1), 3)
F2  = at(sp.diff(F,u,2), 3)
print("F_s =", F_s)
print("F1  =", F1)
print("F2  =", F2)

sig2 = -1/phi2   # ~ eps^{-1}
# gaussian-moment correction terms, each relative to F_s:
t_F2  = sp.simplify(F2*sig2/(2))                 # x^2 moment * F2/2
t_p4  = sp.simplify(phi4*sig2**2*sp.Rational(3,1)/24*F_s)
t_p33 = sp.simplify(phi3**2*sig2**3*sp.Rational(15,1)/72*F_s)
t_F1p3= sp.simplify(F1*phi3*sig2**2*sp.Rational(3,1)/6)
corr  = sp.series(sp.simplify((t_F2 + t_p4 + t_p33 + t_F1p3)/F_s), eps, 0, 3).removeO().expand()
print("relative loop correction (from moments) =", corr)

# also expand F_s itself (contains P(us) and Q(us) eps-terms)
Fs_exp = sp.series(F_s, eps, 0, 3).removeO().expand()
print("F_s expanded =", Fs_exp)
# leading: z_s^{1/4} eps^{-1/2} (1 + (5/4) eps sqrt(zs) - (3/4) eps sqrt(zs) + ...) = ... (1/2) eps sqrt(zs)

# e^{phi_s} correction
print("phi_s eps term coefficient: -zs**(3/2)/12 confirmed:", sp.expand(phi_s).coeff(eps,1))
