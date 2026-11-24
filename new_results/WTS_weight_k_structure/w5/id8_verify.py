# Identification step 8: (1) c_3 closed form vs the D=3300 value of c_3; (2) c_0..c_3 (+delta) reproduce W_5(k) k=4..11 from values_w5.json.
import sys, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 3350                       # set AFTER imports
print("mp.mp.dps =", mp.mp.dps)
R = json.load(open("results/cm_D3300_m11.json"))
c = [mp.mpf(x) for x in R["c"]]; delta = mp.mpf(R["delta"])
lam = mp.log(2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5)
c3cf = 5*z5/64 + 5*pi**2*z3/576 + pi**4*lam/384 + 5*z3*lam**2/96 + 5*pi**2*lam**3/864 + lam**5/288
print("c_3 closed form - c_3(D=3300 data) =", mp.nstr(c3cf - c[3], 5), "(data has 3300 digits)")
V = json.load(open("values_w5.json"))
vals = {int(d["k"]): mp.mpf(d["value"]) for d in V["values"]}
for k in range(4, 12):
    P = sum(c[m]*mp.mpf(k)**m for m in range(4))
    pred = P + (delta if k == 4 else 0)
    print(f"k={k}: |W_5(k) stored - (cubic{' + delta' if k==4 else ''})| = {mp.nstr(abs(vals[k]-pred), 5)}   [stored digits 135]")
