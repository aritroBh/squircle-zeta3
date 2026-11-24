from fractions import Fraction as Fr
import importlib, time
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_ft import ft_delta22

# n=12: theta = 2 sin(pi/12), minpoly x^4-4x^2+1
t0 = time.time()
C12 = Corner("c=sin(pi/12) (n=12)", [1, 0, -4, 0, 1], (Fr(0), Fr(1, 2)), verbose=False)
for lbl, el in [("c", C12.c0), ("1-c", C12.F.sub(C12.one, C12.c0)), ("1+c", C12.F.add(C12.one, C12.c0))]:
    print(f"[{lbl}] = {C12.vec(el)}")
SOm12 = C12.omega_symbol()
d22 = ft_delta22(SOm12, 4, 4)
print(f"n=12 (time {time.time()-t0:.1f}s): TRUE delta_22 #components {len(d22)}")
for k, v in sorted(d22.items())[:6]:
    print("   ", k, v)

# cos(pi/12) = 1/(2 theta) in same field
F = C12.F
c_cos = F.scale(Fr(1, 2), F.inv(F.theta()))
t0 = time.time()
Cc = Corner("c=cos(pi/12) (alpha=pi/6)", [1, 0, -4, 0, 1], c_cos, basis=C12.basis, verbose=False)
for lbl, el in [("c", Cc.c0), ("1-c", Cc.F.sub(Cc.one, Cc.c0)), ("1+c", Cc.F.add(Cc.one, Cc.c0))]:
    print(f"[{lbl}] = {Cc.vec(el)}")
SOmc = Cc.omega_symbol()
d22c = ft_delta22(SOmc, 4, 4)
print(f"cos(pi/12) (time {time.time()-t0:.1f}s): TRUE delta_22 #components {len(d22c)}")
for k, v in sorted(d22c.items())[:6]:
    print("   ", k, v)
