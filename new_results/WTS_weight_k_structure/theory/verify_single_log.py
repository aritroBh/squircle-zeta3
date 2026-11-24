"""Verify (exactly, rational t) the closed form of the generating integrand used in the structural discussion:
   Psi = (1/2) e^{eps*l1} (1/t) [ sqrt(1+y^2) - 1 - y ],   y = w^{-1} e^{eps*(l2-l1)},  w=(1-t)/t,
i.e. the whole M-integrand depends on the two formal logs only through l1 (pure exponential prefactor) and s=l2-l1
(through y).  Compare the coefficient of eps^k l1^m l2^{k-m} with genintegrand's closed form
   c_m(t) = q S_m(t)/(2 t (1-t) q2^k) - delta_{m,k}/(2 k! t) - delta_{m,0}/(2 k! (1-t)),  q = sqrt(t^2+(1-t)^2).
(both are the dt-coefficients, i.e. before multiplication by dt/dv).  k = 1..7, t = 3/10 and t = 7/11 (exact, sympy)."""
import sys, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from math import factorial
from genintegrand import S_poly
eps, l1, l2 = sp.symbols('eps l1 l2')
ok = True
for tval in (sp.Rational(3, 10), sp.Rational(7, 11)):
    t = tval
    w = (1 - t) / t
    s = l2 - l1
    y_expr = sp.exp(eps * s) / w
    Psi = sp.Rational(1, 2) * sp.exp(eps * l1) / t * (sp.sqrt(1 + y_expr**2) - 1 - y_expr)
    K = 7
    ser = sp.series(Psi, eps, 0, K + 1).removeO()
    ser = sp.expand(ser)
    q2 = t**2 + (1 - t)**2
    q = sp.sqrt(q2)
    for k in range(1, K + 1):
        coeff_k = sp.expand(ser.coeff(eps, k))
        for m in range(0, k + 1):
            lhs = sp.simplify(coeff_k.coeff(l1, m).coeff(l2, k - m))
            Sm = S_poly(k, m)
            Sm_val = sum(sp.Rational(c.numerator, c.denominator) * t**e for e, c in enumerate(Sm))
            rhs = q * Sm_val / (2 * t * (1 - t) * q2**k)
            if m == k: rhs -= sp.Rational(1, 2 * factorial(k)) / t
            if m == 0: rhs -= sp.Rational(1, 2 * factorial(k)) / (1 - t)
            d = sp.simplify(lhs - rhs)
            if d != 0:
                ok = False
                print("MISMATCH t=", tval, "k=", k, "m=", m, sp.N(lhs), sp.N(rhs))
print("closed form Psi = (1/2) e^{eps l1}/t [sqrt(1+y^2)-1-y] reproduces genintegrand coefficients for k=1..7, two exact t values:", ok)
