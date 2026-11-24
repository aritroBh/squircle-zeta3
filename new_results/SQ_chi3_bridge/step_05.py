import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

print("Type of POINTS:", type(C.POINTS))
# Let's check the first few points
if hasattr(C.POINTS, 'keys'):
    print("Number of points:", len(C.POINTS))
    for i, k in enumerate(list(C.POINTS.keys())[:10]):
        print(f"  {k}: {C.POINTS[k]}")