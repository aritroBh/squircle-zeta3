# Check what keys the L2mr2h vector uses
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction as Fr

L2_half = C.atom(('li', 2, C.Q(1)/2))
L2r2h = C.atom(('li', 2, C.A / C.Q(2)))
L2mr2h = C.add((L2_half, Fr(1,2)), (L2r2h, -1))

print("L2mr2h vector keys:")
print(list(L2mr2h.keys())[:20])
print(f"\nTotal keys: {len(L2mr2h)}")

# Check if any key contains 4
has_4 = any(4 in w if isinstance(w, tuple) else False for w in L2mr2h.keys())
print(f"Contains key with '4': {has_4}")

# The keys look like they are tuples of integers (Cayley exponents 0-7)
# (4, 0) is NOT in the vals8 dict keys because the vals8 dict has H-word format
# H-words are tuples of length 3, not length 2

# Let's check the format
print(f"\nL2mr2h key types: {[type(k) for k in list(L2mr2h.keys())[:5]]}")
print(f"L2mr2h key lengths: {[len(k) if isinstance(k, tuple) else 'not tuple' for k in list(L2mr2h.keys())[:5]]}")

# Check weight3_with_lifts columns
W = pickle.load(open('research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle', 'rb'))
print(f"\nweight3_with_lifts keys: {list(W.keys())}")
print(f"cols sample: {list(W['cols'])[:10]}")
print(f"Number of cols: {len(W['cols'])}")