import sys, itertools
sys.path.insert(0, ".")
from fractions import Fraction as F
from qs3_symbol import *
import obstruction32 as OB

# candidate S-units z with z and 1-z both in the group
def factors(z):
    try:
        logvec(z); logvec(q3sub(ONE, z)); return True
    except Exception:
        return False
pts = set()
for a in range(-2, 3):
    for b in range(-3, 4):
        for c in range(-2, 3):
            z = q3mul(q3pow(U, a), q3mul(q3pow(V, b), q3pow(R3, c)))
            if z != ZERO and z != ONE and factors(z):
                pts.add(z)
pts = sorted(pts)
print("candidate points:", len(pts))
pairs = []
for x in pts:
    for y in pts:
        try:
            logvec(x); logvec(y); logvec(q3sub(ONE, y)); logvec(q3sub(ONE, x))
            logvec(q3sub(ONE, q3mul(x, y)))
            pairs.append((x, y))
        except Exception:
            pass
print("valid (x,y) pairs:", len(pairs))
pickle = __import__('pickle')
pickle.dump({'pts': pts, 'pairs': pairs}, open("li32_candidates.pkl","wb"))
