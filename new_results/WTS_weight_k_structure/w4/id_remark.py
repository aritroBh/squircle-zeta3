"""Exact (sympy) check of the algebraic rewriting  c_2 = W_2^2 - log2 * S - (7/16) zeta(4),
W_2 = (pi^2+6 log^2 2)/24, S = zeta(3)/8 + pi^2 log2/24 + log^3 2/12, zeta(4) = pi^4/90."""
import sympy as sp
lam, z3, pi = sp.symbols("lam z3 pi")
W2 = (pi**2 + 6*lam**2)/24; S = z3/8 + pi**2*lam/24 + lam**3/12
c2 = -(3*pi**4 + 20*pi**2*lam**2 + 20*lam**4 + 120*z3*lam)/960
print("c2 - (W2^2 - lam*S - 7*pi^4/1440) =", sp.simplify(sp.expand(c2 - (W2**2 - lam*S - sp.Rational(7, 1440)*pi**4))))
