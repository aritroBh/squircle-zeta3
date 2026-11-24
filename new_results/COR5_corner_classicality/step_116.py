import sympy as sp
from sympy import sin, pi, sqrt, Rational, nsimplify
from sympy.polys.numberfields import minimal_polynomial
import mpmath as mp

mp.mp.dps = 80
def embed(alpha, prec=80):
    mp.mp.dps = prec
    x = sp.Symbol('x')
    P = minimal_polynomial(alpha, x)
    coeffs = [mp.mpf(sp.N(c, prec)) for c in sp.Poly(P).all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=200)
    return [mp.log(abs(r)) for r in roots], P

for n in range(3, 13):
    c = sin(pi/n)
    x = sp.Symbol('x')
    P = minimal_polynomial(c, x)
    deg = sp.degree(P)
    cs = sp.nsimplify(c)
    print(f"n={n}: c = {cs}, deg {deg}, minpoly {sp.expand(P)}")
