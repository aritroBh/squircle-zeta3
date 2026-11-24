import sympy as sp
from sympy import sin, pi, sqrt, Rational, nsimplify, AlgebraicNumber, MinimalPolynomial, primorial
import mpmath as mp

mp.mp.dps = 80
def embed(alpha, prec=80):
    """log|sigma(alpha)| for all real embeddings sigma, via mpmath roots of minpoly."""
    mp.mp.dps = prec
    P = MinimalPolynomial(alpha, sp.Symbol('x'))
    coeffs = [mp.mpf(sp.N(c, prec)) for c in sp.Poly(P).all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=200)
    # for totally real fields all roots real; keep all (log modulus)
    return [mp.log(abs(r)) for r in roots], P

def get_field_info(n):
    c = sin(pi/n)
    c = sp.nsimplify(c) if n != 7 and n != 9 and n != 11 else c
    # exact algebraic form
    P = MinimalPolynomial(c, sp.Symbol('x'))
    K_deg = sp.degree(P)
    return c, P, K_deg

for n in range(3, 13):
    c, P, deg = get_field_info(n)
    emb, _ = embed(c)
    print(f"n={n}: c = {sp.nsimplify(c)}, deg {deg}, minpoly {P}")
