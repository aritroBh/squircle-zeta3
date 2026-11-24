# 27_lll_Za_M2.py -- is Za (and V4) a combination of the squircle coefficient M2 and simpler constants?
import mpmath as mp, json
from lll_util import lll_relation, describe
from squircle_consts import consts
mp.mp.dps = 420
d = json.load(open("Q1_sub_400.json"))
M0, M1, M2 = consts()
r2 = mp.sqrt(2); lam, L = mp.log(2), mp.log(1+r2); pi = mp.pi; b = r2-1
Li2 = lambda x: mp.polylog(2, x)
B = [("M2",M2),("1",1),("lam",lam),("L",L),("pi2",pi**2),("lam2",lam**2),("lamL",lam*L),("L2",L**2),("Li2b",Li2(b)),
     ("z3",mp.zeta(3)),("pi2lam",pi**2*lam),("pi2L",pi**2*L),("lam3",lam**3),("lam2L",lam**2*L),("lamL2",lam*L**2),("L3",L**3),
     ("lamLi2b",lam*Li2(b)),("LLi2b",L*Li2(b))]
names, vals = [], []
for n, v in B:
    names += [n]; vals += [v]
    if n != "M2": names += ["r2*"+n]; vals += [r2*v]
for key in ("Za", "V4"):
    t = mp.mpf(d[key])
    rel, res, h = lll_relation([t] + vals, 390)
    expect = 390/(len(vals)+1)
    print("%s: height 10^%.2f (random ~10^%.1f)" % (key, float(mp.log10(h)), expect))
    if mp.log10(h) < 0.6*expect: print("   ", key, "=", describe(rel, names))
