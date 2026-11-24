exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
# G at a1=1, a2=1: univariate in a4
g11 = sp.Poly(G.subs({a1:1, a2:1}), a4)
print("G(1,1,a4) degree:", g11.degree())
roots = sp.polys.polyroots  # use real root isolation
from sympy.polys.rootisolation import dup_isolate_real_roots
coeffs = [int(c) for c in g11.all_coeffs()]
rr = dup_isolate_real_roots(coeffs, sp.ZZ)
print("real roots isolating intervals:")
pos = []
for iv in rr:
    lo, hi = iv
    print("  ", iv, " mid ~", sp.N((lo+hi)/2, 6))
    if lo > 0: pos.append(iv)
print("positive roots:", len(pos))
