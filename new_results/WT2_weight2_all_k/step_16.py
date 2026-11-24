import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from mk_integrands import v_integrand, v, S2

for k in [2,4,5]:
    Rk = v_integrand(k)
    print(f"--- k={k}")
    for key in sorted(Rk, reverse=True):
        R = Rk[key]
        num, den = sp.fraction(R)
        print(key, "den:", sp.factor(den, extension=S2), " numdeg:", sp.Poly(num, v).degree() if num.has(v) else 0)
