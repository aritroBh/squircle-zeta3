# Check how atom works for Li_n
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Check the atom implementation
print("Checking atom implementation for li 3:")
z = C.Q(1) - C.A  # Example point
vec = C.atom(('li', 3, z))
print(f"L3 at {z}: {len(vec)} words")

# The key: atom calls dlog with cross-ratio arguments
# Let's trace through: for ('li', n, z), it calls find_cross(z) to get the frame
# Then builds forms from dlog calls with the frame points

# Let me understand by checking a Li_2
print("\nChecking Li_2 at 1-sqrt(2):")
z = C.Q(1) - C.A
L2 = C.atom(('li', 2, z))
print(f"  {len(L2)} words")

# Now the question: how do we build depth-2 G(a1, a2; 1) where a1, a2 are not cross-ratios?
# The letters 1/x and 1/(xy) need to be evaluated

# For Lam: x = 1 - sqrt(2) = 1 - A*sqrt(2) + 1 = ... 
# Actually: 1 - sqrt(2) where sqrt(2) is not directly in the field
# But the cross ratio values are in Q(A) where A = 1 + sqrt(2)

# x1 = 1 - sqrt(2) = 1 - (A - 1) = 2 - A
print(f"\nx1 = 1 - sqrt(2) as Q(A): = {C.Q(1) - C.A}")

# 1/x1 = 1/(1 - sqrt(2)) = (1 + sqrt(2))/(1 - 2) = -(1 + sqrt(2)) = -A
print(f"1/x1 = -A = {C.Q(-1) * C.A}")

# xy for Lam: x * y = (1 - sqrt(2)) * (-sqrt(2)/2) 
# = -sqrt(2)/2 + 1/2*2 = -sqrt(2)/2 + 1 = 1 - sqrt(2)/2
# So 1/(xy) = 1/(1 - sqrt(2)/2) = 2/(2 - sqrt(2))
# = 2 * (2 + sqrt(2)) / (4 - 2) = 2 + sqrt(2) = A + 1

print(f"\nLam: x = 1 - sqrt(2) = {C.Q(1) - C.A}")
print(f"     y = -sqrt(2)/2 = ... need to express in Q(A)")
print(f"     1/x = -A")
print(f"     1/(xy) = 2 + sqrt(2) = A + 1")

# But is -sqrt(2)/2 expressible in Q(A)?
# -sqrt(2)/2 = -(A - 1)/2 = (1 - A)/2
print(f"\n-sqrt(2)/2 = (1 - A)/2 = {C.Q(1)/2 - C.A/2}")