import sys, itertools, json, time
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
from cor6_sym import sshuffle, rref

# weight-5 quotient T^5/D^5 over 2 letters (free Lie_5 coords), NF coefficients
B5 = list(itertools.product(range(2), repeat=5))
def d5_rows():
    rows = []
    for lu in range(1, 5):
        for u in itertools.product(range(2), repeat=lu):
            for v in itertools.product(range(2), repeat=5 - lu):
                rows.append(sshuffle({u: Fr(1)}, {v: Fr(1)}))
    return rows
t0 = time.time()
D5 = d5_rows()
Dmat5, piv5 = rref([[s.get(w, Fr(0)) for w in B5] for s in D5])
free5 = [j for j in range(len(B5)) if j not in piv5]
print(f"T^5/D^5: dim {len(free5)} (expect 6), free words {[B5[j] for j in free5]} ({time.time()-t0:.0f}s)")

def quotient5(SE):
    """SE {index-word -> NF} -> tuple of 6 NF coords."""
    v = [tuple(SE.get(w, (Fr(0), Fr(0)))) for w in B5]
    for i, pj in enumerate(piv5):
        if any(v[pj]):
            q = v[pj]
            v = [tuple(vv[t] - q[t] * Dmat5[i][j] for t in range(2)) for j, vv in enumerate(v)]
    for j in range(len(B5)):
        if j not in free5:
            assert not any(v[j]), "not integrable"
    return tuple(v[j] for j in free5)
import pickle
pickle.dump({'B5': B5, 'piv5': piv5, 'free5': free5, 'Dmat5': [[str(x) for x in row] for row in Dmat5]},
            open(GOLD + 't5d5.pkl', 'wb'))
print("saved t5d5.pkl")
