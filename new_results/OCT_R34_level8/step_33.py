import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Check the class name
print("Available in cmzv:")
for attr in dir(C):
    if 'S2' in attr or 'qs' in attr.lower() or attr.startswith('Q'):
        print(f"  {attr}")

# Check if a2 = 2+sqrt(2) is in cross_ratios
a2_correct = C.Q(1) + C.A  # This gives (2+1*sqrt2)

print(f"\na2 = 2+sqrt(2) = A+1 = {a2_correct}")

try:
    frame_a2 = C.find_cross(a2_correct)
    print(f"Frame for 2+sqrt(2): {frame_a2}")
except ValueError as e:
    print(f"find_cross(2+sqrt(2)): {e}")