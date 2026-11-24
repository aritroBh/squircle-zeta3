import sympy as sp
# Symbolic verification of the algebra steps (PROVED identities)
a, t = sp.symbols('a t', positive=True)

# (iv) derivative identity: d/da [a^3/sqrt(1+a^4)] = 3a^2(1+a^4)^{-1/2} - 2a^6(1+a^4)^{-3/2}... check:
expr = sp.simplify(sp.diff(a**3/sp.sqrt(1+a**4), a) - 2*a**2/(1+a**4)**sp.Rational(3,2) - a**2/sp.sqrt(1+a**4))
print("(iv) diff(a^3/sqrt(1+a^4)) - 2a^2/(1+a^4)^{3/2} - a^2/sqrt(1+a^4) =", sp.simplify(expr))

# (iii) substitution check: t = sqrt(2) a / sqrt(1+a^4); show  t^2/sqrt(1-t^4) dt = 2 sqrt(2) a^2/(1+a^4)^{3/2} da
tt = sp.sqrt(2)*a/sp.sqrt(1+a**4)
dt_da = sp.diff(tt, a)
integrand_t = tt**2/sp.sqrt(1-tt**4) * dt_da
print("(iii) substituted integrand =", sp.simplify(integrand_t))
print("      target 2sqrt2 a^2/(1+a^4)^{3/2} =", sp.simplify(2*sp.sqrt(2)*a**2/(1+a**4)**sp.Rational(3,2)))
print("      equal:", sp.simplify(integrand_t - 2*sp.sqrt(2)*a**2/(1+a**4)**sp.Rational(3,2)) == 0)

# (ii) Beta evaluation
from sympy import beta, gamma, sqrt, pi
print("(ii) (1/4)B(3/4,1/2) =", sp.simplify(beta(sp.Rational(3,4), sp.Rational(1,2))/4))
print("     = sqrt2 pi^{3/2}/Gamma(1/4)^2 ?", sp.simplify(beta(sp.Rational(3,4),sp.Rational(1,2))/4 - sqrt(2)*pi**sp.Rational(3,2)/gamma(sp.Rational(1,4))**2)==0)
