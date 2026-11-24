# Great! Now I can see: 2+sqrt(2) ≡ (2+1*sqrt2) is in the cross ratios
# And -A ≡ (-1-1*sqrt2) is also there

import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Find the cross for these values
a1_Lam = -C.A  # -A = (-1-1*sqrt2)
a2_Lam_value = C.Q(2) + C.A  # 2 + sqrt(2) = (2+1*sqrt2) 

print(f"Looking for cross-ratio frames:")
print(f"a1 = -A = {a1_Lam}")
print(f"a2 = 2+sqrt(2) = {a2_Lam_value}")

# Try find_cross
try:
    frame1 = C.find_cross(a1_Lam)
    print(f"Frame for -A: {frame1}")
except ValueError as e:
    print(f"find_cross(-A): {e}")

try:
    frame2 = C.find_cross(a2_Lam_value)
    print(f"Frame for 2+sqrt(2): {frame2}")
except ValueError as e:
    print(f"find_cross(2+sqrt(2)): {e}")