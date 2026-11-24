# s11: what "zeta(3)-free" can and cannot mean.  Hexagon with kinks (0, pi/2, pi/2+eps), weights 1:
# as eps -> 0 two kink lines merge and K_p -> the rhombus family {|y|^p + 2|x|^p <= 1}.
# P'''(1) is continuous, the edge-law coefficient is NOT: the jump is carried by the non-local functions,
# whose expansion variable x = r-1 tends to the boundary of convergence (r -> 0), where their sum develops zeta(3).
import mpmath as mp, sympy as sp
from edge_num import P_coeffs_t, zeta3_edge_law, edges
from s04_exact_local import closed_form
from s03_beta_table import z2, z3
mp.mp.dps = 30
ns = [(0, 1), (1, 0)]; dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
out, ok, Zr = closed_form(ns, [1, 2], [0, 0], dirs, "limit rhombus |y|^p + 2|x|^p <= 1")
P3r = mp.mpf(str(sp.N(out[3].subs({z2: sp.pi**2/6, z3: sp.zeta(3)}), 40)))
Zr = mp.mpf(str(sp.N(Zr, 40)))
print(f"limit rhombus: P'''(1) = {mp.nstr(P3r, 25)}, edge law Z = {mp.nstr(Zr, 25)}")
for eps in ('0.1', '0.03', '0.01', '0.003', '0.001'):
    eps = mp.mpf(eps)
    bs = [mp.mpf(0), mp.pi/2, mp.pi/2 + eps]
    P3 = 6*P_coeffs_t(bs, [1, 1, 1], [0, 0, 0])[3]
    Z, parts = zeta3_edge_law(bs, [1, 1, 1])
    rs = []
    for e in edges(bs, [1, 1, 1], [0, 0, 0]):
        i2 = [i for i in range(3) if i not in (e['j0'], e['j1'])][0]
        rs.append(e['c1'][i2]/e['c0'][i2])
    print(f"eps={mp.nstr(eps, 3)}: P'''(1) = {mp.nstr(P3, 20)} (minus limit: {mp.nstr(P3 - P3r, 3)}),"
          f"  Z_hex = {mp.nstr(Z, 20)},  Z_hex - Z_rhombus = {mp.nstr(Z - Zr, 12)},  ratios r = {[mp.nstr(r, 4) for r in rs[:3]]}")
