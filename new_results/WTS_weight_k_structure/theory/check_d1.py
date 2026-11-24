"""d=1 against WT2: the degree-1 generating function is eps*R10 + eps^2/(1-eps)*Rw2.  Exact comparison with genfun.json (Newton data)."""
import json, os
from fractions import Fraction as Fr
J = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "genfun.json")))["1|1,0"]
# R10 = (1-1/sqrt2)/v - 2/(v+a)  : residue 1 - sqrt2/2 at 0, -2 at -a
R10 = {"0": (Fr(1), Fr(-1, 2)), "-a": (Fr(-2), Fr(0))}
# Rw2 = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b)) : residues 1/2 at 0, 1/2 at a, -1 at -a, -1/2 at -b
Rw2 = {"0": (Fr(1, 2), Fr(0)), "a": (Fr(1, 2), Fr(0)), "-a": (Fr(-1), Fr(0)), "-b": (Fr(-1, 2), Fr(0))}
ok = True
for pole in ("0", "a", "-a", "-b"):
    c1 = (Fr(J[pole]["c_d"][0]), Fr(J[pole]["c_d"][1])); D1 = (Fr(J[pole]["D"][0][0]), Fr(J[pole]["D"][0][1]))
    ok &= (c1 == R10.get(pole, (Fr(0), Fr(0)))) and (D1 == Rw2.get(pole, (Fr(0), Fr(0))))
print("d=1: c(1) = residues of R10 and D_1 = residues of Rw2 at all four poles (exact):", ok)
