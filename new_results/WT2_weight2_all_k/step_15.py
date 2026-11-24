import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from mk_integrands import v_integrand, v, S2

t0=time.time()
Rk = v_integrand(3)
print("v_integrand(3) time", time.time()-t0)
for key in sorted(Rk, reverse=True):
    R = Rk[key]
    num, den = sp.fraction(R)
    print(key, "den factors:", sp.factor(den, extension=S2))
