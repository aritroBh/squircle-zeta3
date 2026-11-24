import sys
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L
from fractions import Fraction as Fr

mp.mp.dps = 1120
F4 = mp.mpf(open("c13_W4_hpl_1150.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2
B4_NAMES = ["Li4(-2)", "Li4(1/3)", "Li4(2/3)", "Li4(1/2)", "Li4(1/9)",
            "Li3(1/3)*log2", "Li3(1/3)*log3", "Li3(-1/2)*log2", "Li3(-1/2)*log3",
            "zeta3*log2", "zeta3*log3",
            "Li2(1/3)^2", "Li2(1/3)*log2^2", "Li2(1/3)*log2*log3", "Li2(1/3)*log3^2", "pi2*Li2(1/3)",
            "pi2*log2^2", "pi2*log2*log3", "pi2*log3^2", "pi4",
            "log2^4", "log2^3*log3", "log2^2*log3^2", "log2*log3^3", "log3^4"]
B = [Li(4, mp.mpf(-2)), Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
     l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
     l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
     pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
     l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
cands = {}
for dps in (700, 1000):
    rows = L.find_relation_rows([F4] + B, dps)
    cands[dps] = L.normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(c*v for c, v in zip(cands[dps], [F4] + B))
    print(f"dps {dps}: height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
print("identical:", cands[700] == cands[1000])
a = cands[1000]; a0 = a[0]
print("\nF4 =")
for i in range(1, len(a)):
    if a[i]:
        print(f"  {str(Fr(-a[i], a0)):>10} * {B4_NAMES[i-1]}")
import json
json.dump({"a0": a0, "coeffs": a[1:], "basis": B4_NAMES}, open("c13_F4_relation.json", "w"))
