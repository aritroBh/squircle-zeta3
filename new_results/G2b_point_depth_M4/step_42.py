
import sys, json, re
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/G2_M4alpha")
import os
os.chdir(BASE + "/G2_M4alpha")
from step_10 import load_terms
import g2_fit as gf

# enumerate the distinct k=3 gap triples and k=2 retained words for the 146 signed HPL words
triples = set()
retained3 = set()
pairs2 = set()   # 2-letter words whose Lambda^2 is needed (gap for k=2, retained for k=3)
c0 = Fr(1,2)
words = []
for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
    if len(w) != 5: continue
    sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
    words.append((w, e, q, 1))
    words.append((gf.neg(w), e, q, sgn))
print("signed words:", len(words))
for w, e, q, s in words:
    pts = (0,) + w + (c0,)
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i+k]; ret = w[:i] + w[i+k:]
            ei, ci = pts[i], pts[i+k+1]
            if k == 3:
                triples.add((gap, ei, ci))
                pairs2.add((ret, 0, c0))
            else:
                retained3.add(ret)
                pairs2.add((gap, ei, ci))
print("distinct k=3 gap triples:", len(triples))
print("distinct retained 3-words:", len(retained3))
print("distinct 2-words for Lambda^2:", len(pairs2))
# classify the gap triples: letters strictly inside the interval?
def classify(gap, ei, ci):
    lo, hi = min(ei, ci), max(ei, ci)
    inside = [a for a in gap if lo < a < hi]
    end_eq = (gap[0] == ei, gap[-1] == ci)
    return inside, end_eq
bad = [t for t in triples if classify(*t)[0]]
print("triples with letters strictly inside (PV cases):", len(bad))
for t in sorted(bad)[:20]: print("   ", t)
eq_first = [t for t in triples if t[0][0] == t[1]]
print("triples with first letter == ei:", len(eq_first))
eq_last = [t for t in triples if t[0][-1] == t[2]]
print("triples with last letter == ci:", len(eq_last))
