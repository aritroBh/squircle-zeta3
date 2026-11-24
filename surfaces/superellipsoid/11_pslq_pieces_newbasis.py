# 11_pslq_pieces_newbasis.py -- P1, P2a, P2c in the basis with Li2(b), Li2(-b^4) (and Li2(1/3) allowed).
import mpmath as mp, json, itertools
from pslq_util import find_relation
mp.mp.dps = 320
d = json.load(open("pieces_300.json"))
r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1
alpha = mp.atan(1/r2); Li2 = lambda x: mp.polylog(2, x)
logs = [("lam",lam),("l3",l3),("L",L)]
basis = [("pi^2", pi**2), ("pi*alpha", pi*alpha)] + [(x[0]+"*"+y[0], x[1]*y[1]) for x, y in itertools.combinations_with_replacement(logs, 2)] \
        + [("Li2(b)", Li2(b)), ("Li2(-b^4)", Li2(-b**4)), ("Li2(1/3)", Li2(mp.mpf(1)/3))]
names = [n for n, _ in basis]; vals = [v for _, v in basis]
for key, scale in (("P1", 1), ("P2a", r2), ("P2c", r2)):
    t = mp.mpf(d[key])/scale
    if key == "P1":
        nm, vv = [], []
        for n, v in basis: nm += [n, "r2*"+n]; vv += [v, r2*v]
    else:
        nm, vv = names, vals
    out = find_relation(t, vv, nm, 300, maxcoeff=10**7)
    rel, res = out
    c0 = rel[0]
    print("%s%s = (1/%d)[%s]   residual %s" % (key, "" if scale == 1 else "/sqrt2", -c0,
          " + ".join("%d*%s" % (ci, n) for ci, n in zip(rel[1:], nm) if ci), mp.nstr(res, 3)))
