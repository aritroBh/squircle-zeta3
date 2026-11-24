import json, math
from fractions import Fraction as Ff
J = json.load(open('out/d13/lll_capd17_P1800.json'))
rels = J['rels']
print("n rels saved:", len(rels))
# Load atom names
import pickle
A = pickle.load(open('out/d13/atoms5_cap_d17.pkl','rb'))
print("atom names[:5]:", A['atom_names'][:5], "... total", len(A['atom_names']))
print("keepc:", A['keepc'])
for r in rels[:6]:
    a = [int(x) for x in r['a']]
    g = 0
    for x in a: g = math.gcd(g, abs(x))
    ap = [x//g for x in a]
    mxp = max(abs(x) for x in ap)
    print(f"maxcoef {r['maxcoef']} T-coef {a[0]} gcd {g} -> prim T-coef {ap[0]} prim maxcoef {mxp} nterms {sum(1 for x in ap if x)} resid {r['resid']}")
