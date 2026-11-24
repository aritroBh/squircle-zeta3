from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit
importlib.reload(cor6_nf); importlib.reload(cor6_sunit)
from cor6_nf import NF
from cor6_sunit import *

F = NF([-2, 0, 1])  # Q(sqrt2), theta = sqrt2
S0 = [2]
basis, pool = build_basis(F, S0, 2)
print("basis:", basis)
th = F.theta()
c = (Fr(0), Fr(1, 2))
one = F.one()
for label, el in [("c", c), ("1-c", F.sub(one, c)), ("1+c", F.add(one, c)), ("2", F.from_int(2))]:
    ev = exponent_vector(F, el, basis)
    print(f"[{label}] = {ev}")
