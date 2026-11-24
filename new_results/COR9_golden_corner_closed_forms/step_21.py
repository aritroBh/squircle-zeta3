WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
from cor9_lib import h_plus_order, fund_unit
from sympy.solvers.diophantine.diophantine import diop_DN
import math
D = 103
print("h+(412) =", h_plus_order(412))
x, y = fund_unit(D); print("eps =", x, y, "N =", x*x - D*y*y)
for tgt in (-2, 2, -3, 3, -5, 5, -7, 7):
    sols = diop_DN(D, tgt)
    print(f"x^2-103y^2={tgt}: solutions {sols[:2]}{'...' if len(sols)>2 else ''}")
