from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_ft import ft_delta22

for label, cval in [("c=2/5", Fr(2, 5)), ("c=3/7", Fr(3, 7))]:
    C = Corner(label, [-1, 1], (cval,))
    print(label, "basis:", C.basis)
    SOm = C.omega_symbol()
    d22 = ft_delta22(SOm, 1, C.dimV)
    print(f"  TRUE delta_22(Omega({cval})): #components {len(d22)}")
    for k, v in sorted(d22.items())[:16]:
        print("   ", k, v)
    print()
