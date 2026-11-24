# Check if depth-2 integral produces weight-3 H-words
# The forms_path function builds iterated integrals
# Let's test with a simple example

import sys
import cmzv as C
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')

# G(0, a1, a2; 1) would be forms_path([dlog(0, pinf), dlog(a1, pinf), dlog(a2, pinf)], p0, p1)
# Test: what does forms_path([dlog(0), dlog(1), dlog(2)]) produce?

# First check dlog format
dl0 = C.dlog(C.Q(0), None)
print(f"dlog(0, inf): {dl0}")

dl1 = C.dlog(C.Q(1), None)
print(f"dlog(1, inf): {dl1}")

# Check forms_path with simple forms
forms = [C.dlog(C.Q(0), None), C.dlog(C.Q(1), None), C.dlog(C.A, None)]
vec = C.forms_path(forms, C.Q(0), C.Q(1))
print(f"\nforms_path([dlog(0), dlog(1), dlog(A)], 0, 1):")
print(f"  Keys: {list(vec.keys())[:10]}")
print(f"  Key lengths: {[len(k) for k in list(vec.keys())[:5]]}")
print(f"  Total entries: {len(vec)}")