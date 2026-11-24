import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Check C.POINTS content
print("C.POINTS:")
for k, v in C.POINTS.items():
    print(f"  {k}: exponent {v}")

# Let's see what cross ratios we can form
from itertools import permutations
points_list = [None, C.Q(0), C.Q(1), C.Q(-1), C.A, C.B, -C.A, -C.B]

print("\nEnumerating cross ratios from configuration points:")
cross_vals = {}
for e, p0, p1, pinf in permutations(points_list, 4):
    if pinf is not None and p0 is not None and e is not None:
        # Check if pinf is between p0 and e
        # But we don't have simple comparison for QS2 objects
        cr = C.cross(e, p0, p1, pinf)
        cross_vals[cr] = (e, p0, p1, pinf)

print(f"Number of cross ratios formed: {len(cross_vals)}")
print("Cross ratio values:")
for cr, points in sorted(cross_vals.items(), key=lambda x: str(x[0])):
    print(f"  {cr} from {points}")