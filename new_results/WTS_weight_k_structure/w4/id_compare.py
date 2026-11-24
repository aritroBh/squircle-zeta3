"""Compare the LLL relations (results/id_lll_dps*_D*.json) for equality up to overall sign; print the normalised exact coordinates."""
import json, glob, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
files = sorted(glob.glob(os.path.join(HERE, "results", "id_lll_dps*_D*.json")))
norm = {}
for fn in files:
    d = json.load(open(fn))
    for nm, r in d["results"].items():
        rel = r["relation"]; t = rel["target"]
        # express  target_value = sum (-coef/target_coef) * atom
        norm.setdefault(nm, {})["dps%d_D%d" % (d["dps"], d["D"])] ={k: str(Fr(-v, t)) for k, v in rel.items() if k != "target"}
allok = True
for nm in sorted(norm):
    runs = norm[nm]; ref = next(iter(runs.values()))
    same = all(v == ref for v in runs.values())
    allok &= same
    print(nm, "runs:", sorted(runs), "-> identical:", same, " nterms:", len(ref))
json.dump(norm, open(os.path.join(HERE, "results", "id_lll_normalised.json"), "w"), indent=1)
print("ALL IDENTICAL:", allok)
