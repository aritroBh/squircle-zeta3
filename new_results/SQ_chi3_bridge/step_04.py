import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

print("C.A =", C.A)
print("C.B =", C.B)
print("C.SQ2 =", C.SQ2)
print("\nC.Q =", C.Q)
print("\nPOINTS:")
print(C.POINTS[:10] if hasattr(C.POINTS, '__getitem__') else C.POINTS)