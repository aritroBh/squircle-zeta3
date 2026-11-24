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

w_str = '[0, 1, -1]'
a1, a2, a3 = ast.literal_eval(w_str)
Cw = CD[w_str]
print("C_w built. Differentiating...")
Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])
print("E_w terms:", len(Ew.args))
dE = sp.expand(dpoly(Ew))
ok, coeffs = check_weight1_zero(dE)
print("(i) dE_w/dx log-coefficients all zero:", ok)
for k, v in coeffs.items():
    print(f"    coeff[{k}] = {v}")
