# 18_pslq_Q1_squircle.py -- integer-relation searches for the "squircle-world" pieces of Q1
# (Za, Zb, V4) in the weight<=3 basis over Q(sqrt2) used for M2 in main.tex.
import mpmath as mp, json, sys
from pslq_util import find_relation
mp.mp.dps = 420
d = json.load(open("Q1_sub_400.json"))
r2 = mp.sqrt(2); lam, L = mp.log(2), mp.log(1+r2); pi = mp.pi; b = r2-1
Li2 = lambda x: mp.polylog(2, x); Li3 = lambda x: mp.polylog(3, x)
B = [("1",1),("lam",lam),("L",L),("pi2",pi**2),("lam2",lam**2),("lamL",lam*L),("L2",L**2),("Li2b",Li2(b)),
     ("z3",mp.zeta(3)),("Li3b",Li3(b)),("Li3mb",Li3(-b)),("Li3r",Li3(1/r2)),("pi2lam",pi**2*lam),("pi2L",pi**2*L),
     ("lam3",lam**3),("lam2L",lam**2*L),("lamL2",lam*L**2),("L3",L**3),("lamLi2b",lam*Li2(b)),("LLi2b",L*Li2(b))]
names, vals = [], []
for n, v in B: names += [n, "r2*"+n]; vals += [v, r2*v]
for key in sys.argv[1:] or ["Za", "Zb", "V4", "ZaZb"]:
    t = mp.mpf(d["Za"]) + mp.mpf(d["Zb"]) if key == "ZaZb" else mp.mpf(d[key])
    out = find_relation(t, vals, names, 400, maxcoeff=10**9, maxsteps=10**6)
    if out:
        rel, res = out
        c0 = rel[0]
        print("%s = (1/%d)*[%s]   residual %s" % (key, -c0, " + ".join("%d*%s" % (c, n) for c, n in zip(rel[1:], names) if c), mp.nstr(res, 3)))
    else:
        print(key, ": no relation")
