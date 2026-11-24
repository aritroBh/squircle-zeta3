import sympy as sp
from sympy import sin, pi, sqrt, Rational, nsimplify
from sympy.polys import MinimalPolynomial
import mpmath as mp

mp.mp.dps = 80
def embed(alpha, prec=80):
    mp.mp.dps = prec
    x = sp.Symbol('x')
    P = MinimalPolynomial(alpha, x)
    coeffs = [mp.mpf(sp.N(c, prec)) for c in sp.Poly(P).all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=200)
    return [mp.log(abs(r)) for r in roots], P

for n in range(3, 13):
    c = sin(pi/n)
    x = sp.Symbol('x')
    P = MinimalPolynomial(c, x)
    deg = sp.degree(P)
    emb, _ = embed(c)
    cs = sp.nsimplify(c)
    print(f"n={n}: c = {cs}, deg {deg}")
    print(f"     minpoly: {sp.expand(P)}")
