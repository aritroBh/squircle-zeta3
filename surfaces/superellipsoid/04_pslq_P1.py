# 04_pslq_P1.py -- integer relation search for the 2-D piece P1 (Q(sqrt2) world).
import mpmath as mp, json
from pslq_util import find_relation, show
mp.mp.dps = 300
d = json.load(open("pieces_300.json"))
P1 = mp.mpf(d["P1"])
r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1; Li2b = mp.polylog(2, b); pi = mp.pi
G = mp.catalan
base = [("1",1),("lam",lam),("L",L),("pi2",pi**2),("lam2",lam**2),("lamL",lam*L),("L2",L**2),("Li2b",Li2b)]
names, vals = [], []
for n, v in base:
    names += [n, "r2*"+n]; vals += [v, r2*v]
for dps in (150, 280):
    out = find_relation(P1, vals, names, dps)
    if out:
        rel, res = out
        print("dps", dps, "relation:", rel, " residual", mp.nstr(res, 3))
        print("  P1 =", show(rel, names))
    else:
        print("dps", dps, "no relation")
