import sympy as sp, ast
import mpmath as mp
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

mp.mp.dps = 40
fw = sp.lambdify(x, Ew, 'mpmath')
for xv in [mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('0.8')]:
    print(f"E_w({xv}) =", mp.nstr(fw(xv), 5))

# also check reduced form
Gr = reduce_li2(Ew)
fg = sp.lambdify(x, Gr, 'mpmath')
for xv in [mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('0.8')]:
    print(f"reduce_li2(E_w)({xv}) =", mp.nstr(fg(xv), 5))
