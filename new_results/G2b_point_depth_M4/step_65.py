
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
from g2b_dbg import s05
C0 = Fr(1,2)

needed = {}
for w, coef in s05.signed_words():
    pts = (Fr(0),) + tuple(w) + (C0,)
    for i in range(3):
        gap = w[i:i+3]; ret = w[:i] + w[i+3:]
        ei, ci = pts[i], pts[i+4]
        lam = s05.lam2(ret, C0, Fr(0))
        if lam:
            needed[(gap, ei, ci)] = True
print("needed k=3 gap triples (lam != 0):", len(needed))
# also for the two N5 atoms
for (x, y) in [(Fr(8,9), Fr(3,4)), (Fr(2), Fr(1,3))]:
    w = s05.li32_word(x, y)
    pts = (Fr(0),) + tuple(w) + (Fr(1),)
    for i in range(3):
        gap = w[i:i+3]; ret = w[:i] + w[i+3:]
        ei, ci = pts[i], pts[i+4]
        lam = s05.lam2(ret, Fr(1), Fr(0))
        if lam:
            needed[(gap, ei, ci)] = True
print("with atom gaps:", len(needed))
for t in sorted(needed): print("  ", t)
