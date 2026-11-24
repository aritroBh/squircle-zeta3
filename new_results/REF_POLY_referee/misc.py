"""(i) the mixed terms of a regular-hexagon edge carry more than zeta(3);
(ii) symbolic E3 (edge_sym) against the polar jet, edge by edge, on non-fully-local polygons;
(iii) convention dependence: weighted parallelogram."""
import sympy as sp, pickle, mpmath as mp
from edge_sym import t, L, al, be, B1, C1, x0, y0, x1, y1, X0, X1, Xf1, Xf2
from exact_laws_lib import mom, lt, l1
E = pickle.load(open('edge_E_general.pkl', 'rb'))
E3 = sp.Poly(E[3], X0, X1, Xf1, Xf2)
g21 = E3.coeff_monomial(X0**2*X1); g12 = E3.coeff_monomial(X0*X1**2)
# (i) regular hexagon edge: V0 = R(cos 30, sin 30)?  K_1 of n_i=(cos i pi/3, sin i pi/3): vertex on line <n_3,x>=0 i.e. x=0
R = 1/sp.sqrt(3)
V0 = (0, R); V1 = (-R*sp.sin(sp.pi/3), R*sp.cos(sp.pi/3))
rep = {x0: V0[0], y0: V0[1], x1: V1[0], y1: V1[1], al: sp.Rational(1, 2), be: sp.Rational(1, 2), L: R}
mixed = sp.expand((g21.subs(rep))*lt**2*l1 + (g12.subs(rep))*lt*l1**2)
P = sp.Poly(mixed, t, lt, l1)
val = sp.expand(6*sum(c*mom(n, a, b) for (n, a, b), c in P.terms()))
val = sp.expand(val)
print("(i) regular hexagon, one edge: 6*int(mixed terms) =", val)
print("    x 6 edges, zeta(3) part:", sp.simplify(6*val.coeff(sp.zeta(3))), " (C_3 = -9*sqrt(3)/8)")
rest = sp.expand(6*(val - val.coeff(sp.zeta(3))*sp.zeta(3)))
print("    x 6 edges, the rest:", sp.nsimplify(sp.simplify(rest)), "=", sp.N(rest, 20))
