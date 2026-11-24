# s11b: same family at tiny eps: P'''(1) -> limit rhombus value (slowly, ~ eps log^3 eps), Z_hex -> -24 sqrt5/125,
# so the jump of the edge-law coefficient is Z_rh - Z_hex(0+) = 48sqrt5/25 + 24sqrt5/125 = 264 sqrt5/125.
import mpmath as mp, sympy as sp
from edge_num import P_coeffs_t, zeta3_edge_law
mp.mp.dps = 40
P3r = mp.mpf('-2.161877308317571006751704')  # printed by s11 (exact closed form of the limit rhombus)
s5 = mp.sqrt(5)
print("predicted Z_hex(0+) =", mp.nstr(-24*s5/125, 20), "; jump 264 sqrt5/125 =", mp.nstr(264*s5/125, 20))
for eps in ('1e-5', '1e-7', '1e-9', '1e-12'):
    eps = mp.mpf(eps)
    bs = [mp.mpf(0), mp.pi/2, mp.pi/2 + eps]
    P3 = 6*P_coeffs_t(bs, [1, 1, 1], [0, 0, 0])[3]
    Z, parts = zeta3_edge_law(bs, [1, 1, 1])
    print(f"eps={mp.nstr(eps, 3)}: P'''(1) - P'''_rhombus = {mp.nstr(P3 - P3r, 6)}  (eps*log^3 eps = {mp.nstr(eps*mp.log(eps)**3, 4)}),"
          f"  Z_hex = {mp.nstr(Z, 15)},  Z_hex + 24sqrt5/125 = {mp.nstr(Z + 24*s5/125, 4)}")
