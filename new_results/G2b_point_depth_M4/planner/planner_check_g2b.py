# Planner check of G2b's corrected M4(2pi/3): sum the exact flat-angle Taylor series M4 = sum_n mu_{4,n} c^{2n} at c = 1/2
# (G2's exact table, N = 640), and compare with G2b's 1300-digit value and G2's corner_series value.
import json, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 120
d = json.load(open("<project>/research/08_theory/G2_M4alpha/flat_exact_N640_K4.json"))
print(type(d), (list(d.keys())[:6] if isinstance(d, dict) else len(d)))
z = {"z2": mp.zeta(2), "z3": mp.zeta(3), "z4": mp.zeta(4)}
c2 = mp.mpf(1) / 4
M4 = mp.mpf(0)
for n, terms in d["4"].items():
    for mono, q in terms.items():
        a, b, e = (int(t) for t in mono.strip("()").split(","))
        M4 += mp.mpf(Fr(q).numerator) / Fr(q).denominator * z["z2"]**a * z["z3"]**b * z["z4"]**e * c2**int(n)
g2b = mp.mpf("0.0021017551526063123975140574239779012704445906523150088079")
g2 = mp.mpf("0.00210175515260631239751405742397790127044459065231500880790896811193857636736740094076228543")
print("M4(2pi/3) from the exact flat series, N = 640:", mp.nstr(M4, 60))
print("|series - G2b corrected| =", mp.nstr(abs(M4 - g2b), 3), "(G2b value quoted to 58 digits)")
print("|series - G2 corner_series| =", mp.nstr(abs(M4 - g2), 3))
print("|series - the earlier, buggy 3.0778...| =", mp.nstr(abs(M4 - mp.mpf("3.0778066155833419489833813644364368523986262412469630462066")), 5))
