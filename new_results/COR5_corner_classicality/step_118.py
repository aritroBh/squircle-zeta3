import sympy as sp
from sympy import sin, pi
from sympy.polys.numberfields import minimal_polynomial
import mpmath as mp
from flint import fmpz_mat, fmpz

def field_setup(c):
    """Return (minpoly P of c, roots of P at high precision, degree)."""
    x = sp.Symbol('x')
    P = minimal_polynomial(c, x)
    coeffs = [mp.mpf(sp.N(cf, 250)) for cf in sp.Poly(P).all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=500)
    return P, roots, sp.degree(P)

def embed_in_field(el, c, roots, deg):
    """el as sympy expr in c (or independent); return [log|sigma(el)|] for sigma: c -> root."""
    # express el as polynomial in c
    elp = sp.Poly(sp.nsimplify(sp.expand(el)), sp.Symbol('x')) if not el.free_symbols else el
    # el is a sympy expression in terms of c; substitute c = root numerically
    f = sp.lambdify(sp.Symbol('cvar'), el.subs({c: sp.Symbol('cvar')}) if el.has(c) else el, 'mpmath')
    out = []
    for r in roots:
        v = f(r)
        out.append(mp.log(abs(v)) if v != 0 else mp.mpf(-10)**100)
    return out

def log_vecs_consistent(els, c, prec=200):
    mp.mp.dps = prec
    P, roots, deg = field_setup(c)
    return [embed_in_field(el, c, roots, deg) for el in els], P, deg

n = 3
c = sp.sqrt(3)/2
els = [c, 1 - c, 1 + c]
vecs, P, deg = log_vecs_consistent(els, c)
print("field deg:", deg, "minpoly:", P)
for nm, v in zip(["c", "1-c", "1+c"], vecs):
    print(f"  {nm}: {[mp.nstr(x, 8) for x in v]}")
