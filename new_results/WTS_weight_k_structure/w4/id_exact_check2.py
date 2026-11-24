"""EXACT check: combine the individual LLL closed forms of the 8 constants (results/id_each_const_lll.json) with the exact coefficient
table (id_exact_coeffs.json) and compare, coordinate by coordinate in the 58-dim Q-basis, with the LLL relations of c_2, c_1, c_0 and delta."""
import json, os
from fractions import Fraction as F
from id_lib import COEF            # exact coefficients (p + q sqrt2) of each constant for c0, c1, c2, k3
HERE = os.path.dirname(os.path.abspath(__file__))
each = json.load(open(os.path.join(HERE, "results", "id_each_const_lll.json")))
norm = json.load(open(os.path.join(HERE, "results", "id_lll_normalised.json")))
def closed(n):   # constant n = sum (-coef/target) atom
    rel = each[n]["relation"]; t = rel["target"]
    return {k: F(-v, t) for k, v in rel.items() if k != "target"}
def mulsq(d, p, q):   # (p + q sqrt2) * (dict of atom -> coeff; keys with ' *s2' mean sqrt2 multiples)
    out = {}
    def add(k, v):
        if v: out[k] = out.get(k, F(0)) + v
    for k, v in d.items():
        base = k.replace(" *s2", ""); s2 = k.endswith(" *s2")
        # (p + q r2)(v atom)    or (p + q r2)(v r2 atom) = 2 q v atom + p v r2 atom
        if not s2: add(base, p * v); add(base + " *s2", q * v)
        else: add(base + " *s2", p * v); add(base, 2 * q * v)
    return out
def combo(which):
    tot = {}
    for n, c in COEF.items():
        p, q = c[which]
        for k, v in mulsq(closed(n), p, q).items(): tot[k] = tot.get(k, F(0)) + v
    return {k: v for k, v in tot.items() if v}
for which, nm in (("c2", "c2"), ("c1", "c1"), ("c0", "c0")):
    a = combo(which); b = {k: F(v) for k, v in next(iter(norm[nm].values())).items()}
    print(f"{nm}: exact agreement of (sum_n coef_n * LLL closed form of I_n) with the direct LLL relation of {nm}:", a == b)
a = {}
for k, v in combo("c0").items(): a[k] = a.get(k, 0) + v
for j, mult in ((1, 3), (2, 9)):
    for k, v in combo(f"c{j}").items(): a[k] = a.get(k, 0) + mult * v
for k, v in combo("k3").items(): a[k] = a.get(k, 0) - v
a = {k: v for k, v in a.items() if v}
b = {k: F(v) for k, v in next(iter(norm["delta"].values())).items()}
print("delta: exact agreement:", a == b)
