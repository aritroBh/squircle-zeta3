import sympy as sp
a = sp.symbols('a', positive=True)
# (iii) redo: sqrt((1-a^4)^2) = 1-a^4 for 0<a<1
tt = sp.sqrt(2)*a/sp.sqrt(1+a**4)
integrand_t = sp.simplify(tt**2/sp.sqrt(1-tt**4) * sp.diff(tt, a))
integrand_t = integrand_t.subs(sp.sqrt((1-a**4)**2), 1-a**4)
print("(iii) after sqrt((1-a^4)^2)=1-a^4:", sp.simplify(integrand_t))
print("equal to target:", sp.simplify(integrand_t - 2*sp.sqrt(2)*a**2/(1+a**4)**sp.Rational(3,2))==0)

# (ii) manual Beta chain: B(3/4,1/2) = Gamma(3/4) Gamma(1/2) / Gamma(5/4);
# Gamma(5/4) = Gamma(1/4)/4; Gamma(3/4) Gamma(1/4) = pi/sin(pi/4) = pi sqrt(2).
from sympy import gamma, pi, sqrt, Rational, sin
B = gamma(Rational(3,4))*gamma(Rational(1,2))/gamma(Rational(5,4))
B = B.subs(gamma(Rational(5,4)), gamma(Rational(1,4))/4)
B = B.subs(gamma(Rational(3,4)), pi*sqrt(2)/gamma(Rational(1,4)))
print("(ii) (1/4)B(3/4,1/2) =", sp.simplify(B/4))
