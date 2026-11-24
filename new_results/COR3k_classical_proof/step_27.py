import cor3_common
import sympy as sp
print("patched:", cor3_common._dLi2_at(sp.Rational(1,2), 2), " | ", cor3_common._dLi2_at(-1, 1))
from cor3_common import *
# re-check one weight-2 residual is now clean
c0, resid = const_term(-sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12, 3)
print("C_{1,-1} const:", c0, "| residual:", resid)
