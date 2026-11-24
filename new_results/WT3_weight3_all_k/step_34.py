
import os, pickle
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))
from fractions import Fraction as Fr

def drift_analysis(key):
    d = sum(key)
    ks = sorted(k for (k, kk) in rec if kk == key)
    per_order = {}
    for k in ks:
        for m, c in enumerate(rec[(k,key)]['H'][1], start=1):
            per_order.setdefault(m, {})[k] = c
    print(f"== key {key} (d={d}): orders present {sorted(per_order)}")
    for m in sorted(per_order):
        row = per_order[m]
        # check polynomial growth degree of each QSI component
        comps = {}
        for k, c in row.items():
            for idx, val in enumerate([c.p0, c.p1, c.p2, c.p3]):
                comps.setdefault(idx, {})[k] = val
        desc = []
        for idx, vals in sorted(comps.items()):
            ks2 = sorted(vals)
            ys = [vals[k] for k in ks2]
            # smallest r such that (r+1)-th finite difference is zero
            cur = ys[:]
            r = 0
            while len(cur) > 1 and any(x != 0 for x in cur) and r < 8:
                cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
                r += 1
            deg = r-1 if any(x != 0 for x in ys) else -1
            desc.append(f"p{idx}:deg{deg if deg>=0 else '0(z)'}")
        print(f"  order {m}: {desc}")

for key in [(2,0),(1,1),(3,0),(2,1),(4,0),(3,1),(2,2),(5,0),(4,1),(3,2)]:
    drift_analysis(key)
