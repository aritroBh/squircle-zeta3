
import importlib, dec5d_lib
importlib.reload(dec5d_lib)
import dec5d_lib as D

# Hall duality check for d=2 and d=3
for d in (2,3):
    kit = D.make_kit(D.QS5, D.logvec5, tuple(range(d)))
    det3, det2 = D.hall_duality_check(kit)
    L2r, L3r = D.hall_bases(d)
    print(f"d={d}: dim Lie_2={len(L2r)}, dim Lie_3={len(L3r)}, det M3={det3}, det M2={det2}")

# print the d=2 Hall words to compare with M4c's [X,[X,Y]], [Y,[X,Y]]
L2r, L3r = D.hall_bases(2)
print("\nd=2 Lie_2:", L2r)
print("d=2 Lie_3:", L3r)
