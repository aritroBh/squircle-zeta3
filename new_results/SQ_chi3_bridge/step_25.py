import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle

# Rebuild everything again
X = C.atom(('li', 3, C.B))
Y = C.atom(('li', 3, -C.B))
X_minus_Y = C.add((X, 1), (Y, -1))
LHS = {k: v * 24 for k, v in X_minus_Y.items()}

H1 = {(-1, -1, 1): Fraction(1)}
H3 = {(-1, -1, 3): Fraction(1)}
H5 = {(-1, -1, 5): Fraction(1)}
H7 = {(-1, -1, 7): Fraction(1)}
H_sum = C.add((H1, 1), (H7, 1), (H3, -1), (H5, -1))
H_part = {k: v * 16 for k, v in H_sum.items()}

P2 = C.atom(('pi2',))
Lv = C.atom(('log', C.A))
P2_Lv = C.mul(P2, Lv)
Lv2 = C.mul(Lv, Lv)
Lv3 = C.mul(Lv2, Lv)
RHS = C.add((H_part, 1), (P2_Lv, -3), (Lv3, 4))

# DIFF = LHS - RHS
DIFF = C.add((LHS, 1), (RHS, -1))

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

# Check if all words in DIFF are in vals
print("Checking if all words in DIFF are in vals:")
missing = []
for w in DIFF.keys():
    if w not in vals:
        missing.append(w)
        
print(f"Missing words: {len(missing)}")
if missing:
    print("Missing words:", missing)