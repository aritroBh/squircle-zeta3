from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_sym import W4Quotient
from cor6_ft import ft_coords, ft_delta22

# natural basis for Q(sqrt3): pi2 = 1+sqrt3, pi3 = sqrt3, eps = 2+sqrt3
nat = [(Fr(1), Fr(1)), (Fr(0), Fr(1)), (Fr(2), Fr(1))]
C = Corner("c=sqrt3/2 (triangle, n=3)", [-3, 0, 1], (Fr(0), Fr(1, 2)), basis=nat)
for lbl, el in [("c", C.c0), ("1-c", C.F.sub(C.one, C.c0)), ("1+c", C.F.add(C.one, C.c0)), ("2", C.F.from_int(2))]:
    print(f"[{lbl}] = {C.vec(el)}")
SOm = C.omega_symbol()
d22 = ft_delta22(SOm, C.F.d, C.dimV)
print("TRUE delta_22(Omega(sqrt3/2)) in basis {1+sqrt3, sqrt3, 2+sqrt3}:")
for k, v in sorted(d22.items()):
    print("   ", k, v)
print()
print("COR4 naive (function-level): ((0,1),(0,2)): 32sqrt3/3; ((0,2),(1,2)): -8-16sqrt3/3")
