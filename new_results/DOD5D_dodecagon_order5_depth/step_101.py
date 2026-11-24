import sys, pickle, time
sys.path.insert(0, ".")
from fractions import Fraction as F
from qs3_symbol import *
import obstruction32 as OB

pairs = pickle.load(open("li32_candidates.pkl","rb"))['pairs']
L3 = OB.lie3_basis(3); L2 = OB.lie2_basis(3)
def cobracket_vector(S):
    Psi = OB.psi_tensor(S)
    return tuple(OB.pair(Psi, A, B) for (_, A) in L3 for (_, B) in L2)

t0 = time.time()
vecdata = []
seen = set()
for (a, b) in [(4,1),(3,2),(2,3),(1,4)]:
    for (x, y) in pairs:
        w = (ZERO,)*(b-1) + (q3div(ONE, y),) + (ZERO,)*(a-1) + (q3div(ONE, q3mul(x, y)),)
        S = g_symbol(ZERO, w, ONE)
        v = cobracket_vector(S)
        if not any(v): continue
        key = v
        if key in seen: continue
        seen.add(key)
        vecdata.append(((a,b), x, y, v))
print("distinct nonzero cobracket vectors: %d  (%.0fs)" % (len(vecdata), time.time()-t0))
pickle.dump(vecdata, open("mpl_cobrackets.pkl","wb"))
