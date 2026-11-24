import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction

# Build RHS:
# RHS = 16[H(−1,−1,1) + H(−1,−1,7) − H(−1,−1,3) − H(−1,−1,5)] − 3·P2·Lv + 4·Lv^3

# First part: depth-1 H-words
H1 = {(-1, -1, 1): Fraction(1)}
H3 = {(-1, -1, 3): Fraction(1)}
H5 = {(-1, -1, 5): Fraction(1)}
H7 = {(-1, -1, 7): Fraction(1)}

# Sum: H1 + H7 - H3 - H5
H_sum = C.add((H1, 1), (H7, 1), (H3, -1), (H5, -1))

# Multiply by 16
H_part = {k: v * 16 for k, v in H_sum.items()}

print("H_part has", len(H_part), "words")
print("H_part:", H_part)

# P2 and Lv:
P2 = C.atom(('pi2',))
Lv = C.atom(('log', C.A))

# P2 * Lv:
P2_Lv = C.mul(P2, Lv)
print("\nP2 * Lv has", len(P2_Lv), "words")

# Lv^3:
Lv2 = C.mul(Lv, Lv)
Lv3 = C.mul(Lv2, Lv)
print("Lv^3 has", len(Lv3), "words")

# Full RHS:
# RHS = H_part - 3*P2*Lv + 4*Lv^3
RHS = C.add((H_part, 1), (P2_Lv, -3), (Lv3, 4))

print("\nRHS has", len(RHS), "words")
print("First few words of RHS:")
for i, (k, v) in enumerate(list(RHS.items())[:15]):
    print(f"  {k}: {v}")