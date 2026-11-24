from fractions import Fraction as Fr
import importlib, time
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_ft import ft_delta22

# n=8: theta = 2 sin(pi/8) = 2 cos(3pi/8), minpoly x^4-4x^2+2
t0 = time.time()
C8 = Corner("c=sin(pi/8)=cos(3pi/8) (n=8)", [2, 0, -4, 0, 1], (Fr(0), Fr(1, 2)), verbose=False)
for lbl, el in [("c", C8.c0), ("1-c", C8.F.sub(C8.one, C8.c0)), ("1+c", C8.F.add(C8.one, C8.c0))]:
    print(f"[{lbl}] = {C8.vec(el)}")
SOm8 = C8.omega_symbol()
d22 = ft_delta22(SOm8, 4, 4)
print(f"n=8 (time {time.time()-t0:.1f}s): TRUE delta_22 #components {len(d22)}")
for k, v in sorted(d22.items())[:8]:
    print("   ", k, v)
print("basis:", C8.basis)
