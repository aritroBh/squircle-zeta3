import sympy as sp
from cor3_common import *

# test derivative of Li2 atoms -> polylog(1,.) replacement
for nm in ARG_NAMES:
    d = dpoly(sp.polylog(2, ARGS[nm]))
    print(f"d/dx Li2({nm:>14}) =", sp.cancel(d))
