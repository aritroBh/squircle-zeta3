# 18b_lll_Q1_squircle.py -- LLL version of 18_pslq_Q1_squircle.py (faster).
import mpmath as mp, json, sys
from lll_util import lll_relation, describe
mp.mp.dps = 420
d = json.load(open("Q1_sub_400.json"))
r2 = mp.sqrt(2); lam, L = mp.log(2), mp.log(1+r2); pi = mp.pi; b = r2-1
Li2 = lambda x: mp.polylog(2, x); Li3 = lambda x: mp.polylog(3, x)
B = [("1",1),("lam",lam),("L",L),("pi2",pi**2),("lam2",lam**2),("lamL",lam*L),("L2",L**2),("Li2b",Li2(b)),
     ("z3",mp.zeta(3)),("Li3b",Li3(b)),("Li3mb",Li3(-b)),("Li3r",Li3(1/r2)),("pi2lam",pi**2*lam),("pi2L",pi**2*L),
     ("lam3",lam**3),("lam2L",lam**2*L),("lamL2",lam*L**2),("L3",L**3),("lamLi2b",lam*Li2(b)),("LLi2b",L*Li2(b))]
names, vals = [], []
for n, v in B: names += [n, "r2*"+n]; vals += [v, r2*v]
for key in sys.argv[1:] or ["Za", "Zb", "V4"]:
    t = mp.mpf(d[key])
    for D in (300, 380):
        rel, res, h = lll_relation([t] + vals, D)
        print("%s [D=%d] height=%d residual=%s\n    %s = %s" % (key, D, h, mp.nstr(res, 3), key, describe(rel, names)))
