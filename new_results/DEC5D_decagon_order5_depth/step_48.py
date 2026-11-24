
import dec5d_lib as D
from fractions import Fraction as F
from qsqrt2 import QS2
import sys
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
import w4_symbol as W4

# 1) cross-check logvec2 against w4_symbol.logvec on many values
import random
random.seed(5)
bad = 0
for _ in range(300):
    j = random.randint(-6,6); m = random.randint(-8,8)
    z = (D.SQ2**j) * (D.A2**m) * (QS2(-1) if random.random()<0.5 else QS2(1))
    if D.logvec2(z) != W4.logvec(z): bad += 1; print("MISMATCH", z, D.logvec2(z), W4.logvec(z))
print("logvec2 vs w4_symbol.logvec: 300 random S-units, mismatches =", bad)

# 2) letter differences for the decagon, factored exactly
PHI = D.PHI
letters = {'0': D.QS5(0), '1': D.QS5(1), '-phi': -PHI, 'phi^2': PHI**2}
names = list(letters)
print("\nletter differences a_i - a_j factored (j2, k5, m):")
seen = {}
for a in names:
    for b in names:
        if a >= b: continue
        z = letters[a] - letters[b]
        lv = D.logvec5(z)
        seen[(a,b)] = (z, lv)
        print(f"  {a} - ({b}) = {z}  -> log2^{lv[0]} logs5^{lv[1]} logphi^{lv[2]}")
# group generated: which slots appear?
slots_used = set()
for z, lv in seen.values():
    for i, v in enumerate(lv):
        if v: slots_used.add(i)
print("slots used by letter differences:", slots_used, "(0=log2, 1=logsqrt5, 2=logphi)")
