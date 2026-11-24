import time, pickle
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
mp.mp.dps = 1120
r5 = mp.sqrt(5)
d = pickle.load(open(WD+"golden_B4.pkl", "rb"))
B4_NAMES = d["names"]; B4 = [mp.mpf(v) for v in d["vals"]]
nb = len(B4)
vals = B4 + [r5*x for x in B4]
t0 = time.time()
rows = L7.find_relation_rows(vals, 900)
mp.mp.dps = 1120
pure = []
for r in rows:
    g = L7.normalize(r)
    resid = sum(c*v for c, v in zip(g, vals))
    if abs(resid) < mp.mpf('1e-1000'):
        pure.append(g)
print(f"{time.time()-t0:.0f}s; pure-basis relations found: {len(pure)}")
for g in pure:
    terms = []
    for i in range(nb):
        p, q = g[i], g[i+nb]
        if p or q:
            terms.append(f"({p}{'+' if q>=0 else ''}{q}s5)*{B4_NAMES[i]}")
    print(f"  h{max(abs(c) for c in g)}: {' + '.join(terms)[:220]}")
pickle.dump([list(map(int, g)) for g in pure], open(WD+"golden_B4_relations.pkl", "wb"))
