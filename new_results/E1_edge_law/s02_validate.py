# s02: validate the radial-projection integrand E_k(t) against D15b's polar coefficient integrals
# (poly_coeffs.make_coeffs, an independent derivation), edge by edge, and check the closed-form edge law
# against C_m and the rhombus.
import sys, json, mpmath as mp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
import poly_coeffs as pc
from edge_num import *

mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
out = {}
def Cm(m):
    x = mp.pi/(2*m)
    return 48*m*mp.sin(x)**4*mp.cos(x)*(3*mp.sin(x)**2 - 1)

cases = {
  'reg_m2': ([mp.mpf(0), -mp.pi/2], [1, 1], [0, 0]),
  'reg_m3': ([mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 1, 1], [0, 0, 0]),
  'reg_m4': ([-(i*mp.pi/4) for i in range(4)], [1]*4, [0]*4),
  'hexA_w123': ([mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3], [0, 0, 0]),
  'hexB_angles': ([mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1], [0, 0, 0]),
}
for name, (bs, ws, lams) in cases.items():
    E = edges(bs, ws, lams)
    coeffs = pc.make_coeffs(bs, ws, lams)
    worst = [mp.mpf(0)]*4
    totT = [mp.mpf(0)]*4
    for e in E:
        It = [mp.quad(lambda t: E_series(e, t)[k], [0, mp.mpf(1)/2, 1]) for k in range(4)]
        a, b = e['th0'], e['th1']
        if b < a: b += 2*mp.pi
        Ith = [mp.quad(lambda th: coeffs(th)[k], [a, (a+b)/2, b]) for k in range(4)]
        for k in range(4):
            worst[k] = max(worst[k], abs(It[k] - Ith[k]))
            totT[k] += It[k]
    Z, parts = zeta3_edge_law(bs, ws)
    print(name, "edges:", len(E), " per-edge |int_t E_k - int_theta D_k| max, k=0..3:",
          [mp.nstr(x, 3) for x in worst])
    print("   P(1), P'(1), P''(1), P'''(1) =", [mp.nstr(mp.factorial(k)*totT[k], 25) for k in range(4)])
    print("   edge-law zeta(3) coefficient Z =", mp.nstr(Z, 30))
    out[name] = dict(P=[str(mp.factorial(k)*totT[k]) for k in range(4)], Z=str(Z),
                     worst=[str(x) for x in worst])
print("C_2 =", mp.nstr(Cm(2), 30), " 6*sqrt2 =", mp.nstr(6*mp.sqrt(2), 30))
print("C_3 =", mp.nstr(Cm(3), 30), " -9sqrt3/8 =", mp.nstr(-9*mp.sqrt(3)/8, 30))
print("C_4 =", mp.nstr(Cm(4), 30))
# rhombus |x/a|+|y/b| <= 1: normals e2 (kink line y=0, b=0) weight 1/b... kink b=0 is direction x-axis:
for (A, B) in [(1, 2), (3, 1), (mp.mpf('0.7'), mp.mpf('1.9'))]:
    bs = [mp.mpf(0), mp.pi/2]
    # s_j = |sin(th-b_j)|: b=0 -> |y| (normal e2), b=pi/2 -> |x| (normal -e1)
    ws = [1/mp.mpf(B), 1/mp.mpf(A)]
    Z, parts = zeta3_edge_law(bs, ws)
    r = mp.sqrt(mp.mpf(A)**2 + mp.mpf(B)**2)
    print("rhombus a,b =", A, B, " Z =", mp.nstr(Z, 30), " 24a^2b^2/r^3 =", mp.nstr(24*A**2*B**2/r**3, 30),
          " diff", mp.nstr(Z - 24*A**2*B**2/r**3, 3))
json.dump(out, open(f's02_out_{mp.mp.dps}.json', 'w'), indent=1)
