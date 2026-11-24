import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, log, beta
mp.dps = 60
j = 1j
# v(u) = sqrt((1-e^{-2u})/(2u)): Taylor coefficients
import sympy as sp
u = sp.Symbol('u')
v2 = (1-sp.exp(-2*u))/(2*u)
v = sp.sqrt(v2)
ser = sp.series(v, u, 0, 8).removeO()
print("v(u) =", ser)
c = [ser.coeff(u, k) for k in range(6)]
print("c_k =", c)
# g-coefficients: g_k = c_k * Gamma(k+3/2)/Gamma(3/2)
for k in range(1,4):
    g_k = sp.Rational(sp.str(sp.factorial(k))) and c[k]*sp.gamma(sp.Rational(k)+sp.Rational(3,2))/sp.gamma(sp.Rational(3,2))
    print(f"g_{k} =", sp.simplify(c[k]*sp.gamma(sp.Rational(k)+sp.Rational(3,2))/sp.gamma(sp.Rational(3,2))), "=", sp.N(c[k]*sp.gamma(sp.Rational(k)+sp.Rational(3,2))/sp.gamma(sp.Rational(3,2)),10))
