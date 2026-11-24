import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let me understand the POINTS mapping better
print("POINTS dict:")
for k, v in C.POINTS.items():
    print(f"  {k}: {v}")
    
# The POINTS dict maps geometric points to indices
# None -> 0, 0 -> 4, 1 -> 6, -1 -> 2

# Let's check what C.A and C.B are again
print("\nC.A =", C.A)
print("C.B =", C.B)

# Let's check what cross ratios these correspond to
# C.B = √2 - 1 ≈ 0.414...
# C.A = 1 + √2 ≈ 2.414...