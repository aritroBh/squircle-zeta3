from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_ft import ft_delta22

# n=10: c = sin(pi/10) = (sqrt5-1)/4 in Q(sqrt5)
C10 = Corner("c=sin(pi/10) (n=10)", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4)))
print("basis:", C10.basis)
for lbl, el in [("c", C10.c0), ("1-c", C10.F.sub(C10.one, C10.c0)), ("1+c", C10.F.add(C10.one, C10.c0))]:
    print(f"[{lbl}] = {C10.vec(el)}")
SOm10 = C10.omega_symbol()
d22 = ft_delta22(SOm10, C10.F.d, C10.dimV)
print("TRUE delta_22(Omega(sin pi/10)): #components:", len(d22))
for k, v in sorted(d22.items())[:10]:
    print("   ", k, v)
