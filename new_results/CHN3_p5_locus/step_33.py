exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
g11 = sp.Poly(G.subs({a1:1, a2:1}), a4)
from sympy.polys.rootisolation import dup_isolate_real_roots
coeffs = [int(c) for c in g11.all_coeffs()]
rr = dup_isolate_real_roots(coeffs, sp.ZZ)
print("raw:", rr)
# refine numerically
import mpmath as mp
mp.mp.dps = 40
g11f = sp.lambdify(a4, G.subs({a1:1,a2:1}), 'mpmath')
# find roots via polyroots
cf = [mp.mpf(int(c)) for c in g11.all_coeffs()]
roots = mp.polyroots(cf)
print("\nnumeric roots of G(1,1,a4):")
for r in roots:
    if abs(mp.im(r)) < mp.mpf('1e-30') and mp.re(r) > 0:
        print("  positive real:", mp.nstr(mp.re(r), 30))
