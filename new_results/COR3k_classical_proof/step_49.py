import sympy as sp, ast, time
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}
print("=== [0, 1, -1]: E_w == 0 after reduction ===")
t0 = time.time()
a1, a2, a3 = 0, 1, -1
Ew = sp.expand(dpoly(CD['[0, 1, -1]']) - f_a(a1)*C2[(a2, a3)])
ok, terms = check_weight2_zero(Ew)
print("all coefficients 0:", ok, f"({time.time()-t0:.1f}s), {len(terms)} monomials")
if not ok:
    for m, c in terms:
        if c != 0: print("   NONZERO:", m, "=", c)
