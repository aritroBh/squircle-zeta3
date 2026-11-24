# 06_pslq_P2.py -- integer relation searches for P2a/sqrt2 and P2c/sqrt2.
import mpmath as mp, json, itertools
from pslq_util import find_relation, show
mp.mp.dps = 320
d = json.load(open("pieces_300.json"))
P2a = mp.mpf(d["P2a"]); P2c = mp.mpf(d["P2c"])
r2, r3 = mp.sqrt(2), mp.sqrt(3)
lam, l3, L, K, M = mp.log(2), mp.log(3), mp.log(1+r2), mp.log(2+r3), mp.log(r2+r3)
pi = mp.pi
alpha = mp.atan(1/r2)
Li2 = lambda x: mp.polylog(2, x)
logs = [("lam",lam),("l3",l3),("L",L),("K",K)]
prods = [("pi2", pi**2)] + [(a[0]+"*"+b[0], a[1]*b[1]) for a, b in itertools.combinations_with_replacement(logs, 2)]
li2s = [("Li2(1/3)",Li2(mp.mpf(1)/3)),("Li2(b)",Li2(r2-1)),("Li2(r3/2)",Li2(r3/2)),
        ("Li2(4r3-6)",Li2(4*r3-6)),("Li2(2r2/3)",Li2(2*r2/3))]
basis = prods + li2s
names = [n for n, _ in basis]; vals = [v for _, v in basis]
for tname, t in (("P2a/sqrt2", P2a/r2),):
    for dps in (200, 300):
        out = find_relation(t, vals, names, dps)
        print(tname, "dps", dps, ":", (out[0], mp.nstr(out[1], 3)) if out else None)
        if out: print("   ", tname, "=", show(out[0], names))
# P2c: allow pi*alpha
basis_c = prods + [("pi*alpha", pi*alpha), ("alpha^2", alpha**2)] + li2s
names_c = [n for n, _ in basis_c]; vals_c = [v for _, v in basis_c]
for dps in (200, 300):
    out = find_relation(P2c/r2, vals_c, names_c, dps)
    print("P2c/sqrt2 dps", dps, ":", (out[0], mp.nstr(out[1], 3)) if out else None)
    if out: print("    P2c/sqrt2 =", show(out[0], names_c))
