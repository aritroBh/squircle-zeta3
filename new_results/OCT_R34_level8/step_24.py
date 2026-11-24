# Check weight3 pickle cols
import pickle
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')

W = pickle.load(open('research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle', 'rb'))
print(f"weight3_with_lifts keys: {list(W.keys())}")
print(f"cols sample: {list(W['cols'])[:10]}")
print(f"Number of cols: {len(W['cols'])}")

# The values8 pickle has 576 entries, which is different from Li_2 vectors
# Need to check if the depth-2 G-integrals produce weight-3 H-words (length 3)
# Let's check if path produces length-3 tuples

import cmzv as C
# Check a weight-3 atom
L3x1 = C.atom(('li', 3, C.Q(1) - C.A))
print(f"\nL3x1 key lengths: {[len(k) for k in list(L3x1.keys())[:5]]}")
print(f"L3x1 sample keys: {list(L3x1.keys())[:10]}")