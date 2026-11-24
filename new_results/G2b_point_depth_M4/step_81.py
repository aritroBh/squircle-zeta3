
import json, mpmath as mp
mp.mp.dps = 200
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
G2 = "<project>/research/08_theory/G2_M4alpha/"
V = json.load(open(G2B + "values_hex.json"))
M4_hex = mp.mpf(V['M4_hex'])
qline = [l for l in open(G2 + "values4_2pi_3.txt") if l.startswith("M4 =")][0]
M4_quad = mp.mpf(qline.split("=")[1].strip())
print("M4(2pi/3) B4-series (1300 dig):", mp.nstr(M4_hex, 40))
print("M4(2pi/3) G2 quadrature       :", mp.nstr(M4_quad, 40))
print("difference:", mp.nstr(abs(M4_hex - M4_quad), 5))
with open(G2B + "progress.log", "a") as f:
    f.write(f"M4(2pi/3) cross-check vs G2 quadrature: diff = {mp.nstr(abs(M4_hex - M4_quad), 5)}\n")
