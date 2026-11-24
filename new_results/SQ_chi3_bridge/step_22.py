import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Regularize the DIFF vector
DIFF_reg = C.regularize(DIFF)

print("DIFF_reg has", len(DIFF_reg), "words")
print("\nAll words in DIFF_reg:")
for k, v in sorted(DIFF_reg.items()):
    print(f"  {k}: {v}")