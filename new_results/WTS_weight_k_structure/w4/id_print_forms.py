"""Print the normalised exact closed forms (coefficient of each atom as p + q sqrt2) for c_2, c_1, c_0."""
import json, os
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
norm = json.load(open(os.path.join(HERE, "results", "id_lll_normalised.json")))
for nm in ("c2", "c1", "c0"):
    rel = {k: F(v) for k, v in next(iter(norm[nm].values())).items()}
    atoms = sorted({k.replace(" *s2", "") for k in rel})
    print(f"\n{nm} =")
    for a in atoms:
        p = rel.get(a, F(0)); q = rel.get(a + " *s2", F(0))
        print(f"   [{a}]  {p}" + (f"  + ({q}) sqrt2" if q else ""))
