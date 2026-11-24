import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# 2 + sqrt(2) = A + 1 where A = 1 + sqrt(2)
a2_Lam_expected = C.A + C.Q(1)
print(f"2 + sqrt(2) = A + 1 = {a2_Lam_expected}")

# Check if (2+1*sqrt2) is in cross ratios
from itertools import permutations
points_list = [None, C.Q(0), C.Q(1), C.Q(-1), C.A, C.B, -C.A, -C.B]

print("\nSearching for frame giving 2+sqrt(2):")
found = False
for e, p0, p1, pinf in permutations(points_list, 4):
    if pinf is not None and p0 is not None and e is not None and C.Q(0) is not None:
        try:
            cr = C.cross(e, p0, p1, pinf)
            cr_val = repr(cr)
            if "(2+1*sqrt2)" in cr_val or "2+sqrt" in cr_val:
                print(f"  {cr} from e={e}, p0={p0}, p1={p1}, pinf={pinf}")
                found = True
        except:
            pass

# Let me check the enumerated cross values correctly
print("\nComputing cross values again:")
cross_map = {}
for e, p0, p1, pinf in permutations(points_list, 4):
    if pinf is not None and p0 is not None and e is not None and p1 is not None:
        try:
            cr = C.cross(e, p0, p1, pinf)
            cross_map[cr] = (e, p0, p1, pinf)
        except:
            pass

# Check if a2_Lam_expected is in cross_map
if a2_Lam_expected in cross_map:
    print(f"Found frame for {a2_Lam_expected}: {cross_map[a2_Lam_expected]}")
else:
    print(f"{a2_Lam_expected} NOT in cross ratios")
    # Check if it equals any cross ratio
    for cr in cross_map:
        if repr(cr).replace(" ", "") == repr(a2_Lam_expected).replace(" ", ""):
            print(f"  But found matching: {cr}")
            break
    
    # Manual check: compute the numeric value
    import mpmath as mp
    mp.mp.dps = 50
    r2 = mp.sqrt(2)
    a2_mp = 2 + r2
    print(f"\n  Numeric: 2+sqrt(2) = {mp.nstr(a2_mp, 20)}")
    
    # Check each cross ratio value numerically
    print("  Cross ratio values that match numerically:")
    for cr, frame in cross_map.items():
        # cr is a QS2 object, need to evaluate numerically
        # QS2 has .a and .b for rational and sqrt(2) parts
        cr_mp = mp.mpf(cr.a) + mp.mpf(cr.b) * r2
        if abs(cr_mp - a2_mp) < 1e-10:
            print(f"    {cr} = {mp.nstr(cr_mp, 20)} from {frame}")