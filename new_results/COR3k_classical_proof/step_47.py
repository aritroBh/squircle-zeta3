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
Ew = sp.expand(dpoly(CD['[0, 1, -1]']) - f_a(a1)*C2[(a2, a3)])
print("polylogs in Ew:")
for pl in set(t for t in sp.preorder_traversal(Ew) if isinstance(t, sp.polylog)):
    print("   ", pl.args[0], "arg:", pl.args[1], "| key:", arg_key(pl.args[1]))
e0 = reduce_li2(Ew)
print("\npolylogs after reduce_li2:")
for pl in set(t for t in sp.preorder_traversal(e0) if isinstance(t, sp.polylog)):
    print("   ", pl.args[0], "arg:", pl.args[1])
