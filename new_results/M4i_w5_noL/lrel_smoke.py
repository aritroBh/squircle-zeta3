
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path("<project>/research/08_theory/M4i_w5_noL")
sys.path.insert(0, str(HERE))
import mpmath as mp
import m4f_common as MC
import nice_search as NS
import w5_symbol as W
from w5_search import evaluate

digits = 1200
d = MC.build(digits)   # sets prec digits+120, loads cached 65 S0 atom values
vals = d["vals"]; Lchi = d["Lchi"]; ex = d["extras"]
r2 = mp.sqrt(2)
sol = NS.solve(list(NS.S0), "smoke")
names = sol["names"]
recmap = {nm: rec for nm, _, rec in sol["basis"]}
extrasQ = [e["x"] for e in ex]
# atom values: 65 S0 (cached), Li5 extras (fresh polylog at digits+60), Lchi
li5extra = {}
for k, z in enumerate(extrasQ):
    li5extra[f"Li5({z})"] = mp.mpf(evaluate(("li",5,z), digits+60).mid().str(digits+110, radius=False))
def val_of(nm):
    if nm == "Lchi8_5": return Lchi
    if nm in recmap: 
        i = names.index(nm); return vals[i]
    return li5extra[nm]
terms = json.loads((HERE/"L_relation.json").read_text())["terms"]
tot = mp.mpf(0)
for nm,(p,q) in terms.items():
    fp, fq = F(p), F(q)
    tot += (mp.mpf(fp.numerator)/fp.denominator + r2*mp.mpf(fq.numerator)/fq.denominator)*val_of(nm)
res = tot - Lchi
print("smoke test @%d digits:  RHS - L(chi8,5) = %s" % (digits, mp.nstr(res, 8)))
print("Lchi = %s" % mp.nstr(Lchi, 30))
