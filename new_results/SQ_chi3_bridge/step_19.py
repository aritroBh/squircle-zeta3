import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let me build the LHS and RHS vectors
# First, LHS = 24 * (X - Y)

X = C.atom(('li', 3, C.B))
Y = C.atom(('li', 3, -C.B))

# (X - Y):
X_minus_Y = C.add((X, 1), (Y, -1))

# LHS = 24 * (X - Y)
LHS = {k: v * 24 for k, v in X_minus_Y.items()}

print("LHS vector has", len(LHS), "words")
print("First few words of LHS:")
for i, (k, v) in enumerate(list(LHS.items())[:10]):
    print(f"  {k}: {v}")