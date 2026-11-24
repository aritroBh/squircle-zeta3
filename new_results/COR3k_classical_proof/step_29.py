import sympy as sp, ast
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
# find polylog(2,.) terms in E_w
pl2 = [t for t in sp.preorder_traversal(Ew) if isinstance(t, sp.polylog)]
print("polylogs in E_w:", set(pl2))
# terms of dE containing log(x)*log(x+1):
dE = sp.expand(dpoly(Ew))
import sympy
bad = [t for t in dE.args if t.has(sp.log(x)*sp.log(x+1)) or (t.has(sp.log(x)) and t.has(sp.log(x+1)))]
print("\nterms in dE with both log x and log(1+x):")
for t in bad[:5]:
    print("   ", t)
