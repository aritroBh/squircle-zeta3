import sympy as sp
import mpmath as mp
from cor3_common import *

mp.mp.dps = 40
# check each li2_expansion at a small x numerically
xv = sp.Rational(1, 1000)
print("Li2 atom expansions at x=0.001 (error should be ~1e-12 for n=4):")
for nm in ARG_NAMES:
    ser = li2_expansion(nm, 4)
    val_sym = ser.subs({x: xv, LX: sp.log(xv)})
    val = mp.mpf(sp.N(val_sym, 40))
    ref = mp.polylog(2, mp.mpf(str(ARGS[nm].subs(x, xv))))
    print(f"  {nm:>14}: err = {mp.nstr(abs(val-ref), 2)}")

print("\nLi3 atom expansions at x=0.001 (n=2):")
for nm in ARG_NAMES:
    ser = li3_expansion(nm, 2)
    val_sym = ser.subs({x: xv, LX: sp.log(xv)})
    val = mp.mpf(sp.N(val_sym, 40))
    ref = mp.polylog(3, mp.mpf(str(ARGS[nm].subs(x, xv))))
    print(f"  {nm:>14}: err = {mp.nstr(abs(val-ref), 2)}")
