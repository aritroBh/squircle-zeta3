from fractions import Fraction as Fr
import importlib, time
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_ft import ft_delta22

C5 = Corner("c=sin(pi/5) (n=5)", [5, 0, -5, 0, 1], (Fr(0), Fr(1, 2)), verbose=False)
t0 = time.time()
for lbl, el in [("c", C5.c0), ("1-c", C5.F.sub(C5.one, C5.c0)), ("1+c", C5.F.add(C5.one, C5.c0)), ("2", C5.F.from_int(2))]:
    print(f"[{lbl}] = {C5.vec(el)}")
t0 = time.time()
SOm5 = C5.omega_symbol()
print("symbol time", time.time() - t0, "words", len(SOm5))
t0 = time.time()
d22 = ft_delta22(SOm5, 4, 5)
print("delta22 time", time.time() - t0)
print("TRUE delta_22(Omega(sin pi/5)): #components", len(d22))
for k, v in sorted(d22.items()):
    print("   ", k, v)
