# Small exact/numeric side checks used in identify_w5.md:
#  (1) alternative zeta-form of c_3;  (2) c_3 = (5/3) W_2 S + R ;  (3) pole-0 identity  3A22 - 4A31 + A40 = (1/2) int_b^a (dv/v) log^4(t/(1-t)).
import sys, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 60                      # set AFTER imports
print("mp.mp.dps =", mp.mp.dps)
R = json.load(open("results/cm_D3300_m11.json")); c3 = mp.mpf(R["c"][3])
lam = mp.log(2); pi = mp.pi; z2 = mp.zeta(2); z3 = mp.zeta(3); z4 = mp.zeta(4); z5 = mp.zeta(5)
alt = 5*z5/64 + 5*z2*z3/96 + 5*z3*lam**2/96 + 15*z4*lam/64 + 5*z2*lam**3/144 + lam**5/288
print("(1) zeta-form - c3 =", mp.nstr(alt - c3, 5))
W2 = (pi**2 + 6*lam**2)/24; S = z3/8 + pi**2*lam/24 + lam**3/12
Rr = 5*z5/64 - pi**4*lam/3456 - 5*pi**2*lam**3/288 - lam**5/32
print("(2) (5/3) W2 S + R - c3 =", mp.nstr(mp.mpf(5)/3*W2*S + Rr - c3, 5))
B = json.load(open("results/blocks_D3300_m11.json"))["blocks"]
P0 = 3*mp.mpf(B["2,2|0"]["mid"]) - 4*mp.mpf(B["3,1|0"]["mid"]) + mp.mpf(B["4,0|0"]["mid"])
sq2 = mp.sqrt(2); a = 1 + sq2; b = sq2 - 1
f = lambda v: mp.log((v*v + 2*v - 1)/(1 + 2*v - v*v))**4 / v / 2
I = mp.quad(f, [b, 1, a])
print("(3) P0 =", mp.nstr(P0, 30), " (1/2) int dv/v log^4(t/(1-t)) =", mp.nstr(I, 30), " diff =", mp.nstr(P0 - I, 5), "(quadrature limited to ~ working precision)")
