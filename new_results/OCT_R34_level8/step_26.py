# Check forms_path output
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Test forms_path with simple forms
forms = [C.dlog(C.Q(0), None), C.dlog(C.Q(1), None), C.dlog(C.A, None)]
vec = C.forms_path(forms, C.Q(0), C.Q(1))
print(f"forms_path([dlog(0), dlog(1), dlog(A)], 0, 1):")
print(f"  Keys: {list(vec.keys())[:10]}")
print(f"  Key lengths: {set(len(k) for k in vec.keys())}")
print(f"  Total entries: {len(vec)}")

# Now try with the letters for Lam
# Lam letters: a1 = -A, a2 = 2+sqrt(2) 
# But 2+sqrt(2) = A+1 may not be directly a cross-ratio point

# Let me check what POINTS are available
print(f"\nC.POINTS: {C.POINTS}")

# The letters -A and A+1 need to be accessible
# -A is in POINTS (key = 1)
# A+1 = 2+sqrt(2) ... let's check cross ratios