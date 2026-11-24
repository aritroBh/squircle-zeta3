import sympy as sp, ast
import cor3_common
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}
a1, a2, a3 = 0, 1, -1
Cw = CD['[0, 1, -1]']
Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])

PIS = sp.symbols('PIS')
e0 = reduce_li2(Ew)
e0 = e0.subs(dict(zip(LI2_BASIS, cor3_common._ZS)), simultaneous=True)
es = sp.expand(e0).subs(dict(zip(LOG_ATOMS, cor3_common._YS)))
es = es.subs(sp.pi**2, PIS)
gens = list(cor3_common._ZS) + list(cor3_common._YS) + [PIS]
p = sp.Poly(es, *gens, domain=sp.QQ.frac_field(x))
print("total monomials:", len(p.terms()))
for m, c in p.terms():
    cc = sp.cancel(sp.sympify(c))
    if cc != 0:
        print("NONZERO", m, "->", cc)
print("done")
