import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let me check both 2+sqrt(2) = A+1 and the error showed 3+sqrt(2)
print(f"A = {C.A}")
print(f"A+1 = {C.A + C.Q(1)}")
print(f"2+sqrt(2) expected: should be A+1 = 1+sqrt(2)+1 = 2+sqrt(2)")

# Check the QS2 representation
a2_expected = C.A + C.Q(1)
print(f"\na2 = {a2_expected}"  )

# And check that my earlier computation had a bug
a2_wrong = C.Q(2) + C.A
print(f"C.Q(2) + C.A = {a2_wrong}")

# These should be the same. Check:
print(f"Equal? {a2_expected == a2_wrong}")

# Let me trace through: A = 1 + sqrt(2)
# C.A = QS2(rational=1, sqrt2_coeff=1) representing 1 + 1*sqrt(2)
# C.Q(2) = QS2(rational=2, sqrt2_coeff=0)
# C.Q(2) + C.A = QS2(2+1, 0+1) = QS2(3, 1) = 3 + sqrt(2)

# That's wrong! We want 2 + sqrt(2) = QS2(2, 1)
a2_correct = C.QS2(2, 1)
print(f"\nCorrect a2 = QS2(2, 1) = {a2_correct}")

# Now check find_cross
try:
    frame_a2 = C.find_cross(a2_correct)
    print(f"Frame for 2+sqrt(2): {frame_a2}")
except ValueError as e:
    print(f"find_cross(2+sqrt(2)): {e}")