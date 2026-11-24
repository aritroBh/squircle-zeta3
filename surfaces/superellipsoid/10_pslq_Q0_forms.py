# 10_pslq_Q0_forms.py -- closed form of Q0 and s2 in several Li2 bases; pick the simplest.
import mpmath as mp, json, itertools
from pslq_util import find_relation
mp.mp.dps = 320
d = json.load(open("pieces_300.json"))
Q0 = mp.mpf(d["Q0"])
r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1
alpha = mp.atan(1/r2)
Li2 = lambda x: mp.polylog(2, x)
from squircle_consts import consts
M0, M1, M2 = consts()
s2 = 48*M1 - 4*pi**2 - 8*Q0
print("Q0 =", mp.nstr(Q0, 60)); print("s2 =", mp.nstr(s2, 60))
logs = [("lam",lam),("l3",l3),("L",L)]
prods = [("pi^2", pi**2), ("pi*alpha", pi*alpha)] + [(x[0]+"*"+y[0], x[1]*y[1]) for x, y in itertools.combinations_with_replacement(logs, 2)]
A = {"Li2(1/3)": Li2(mp.mpf(1)/3), "Li2(-1/2)": Li2(-mp.mpf(1)/2), "Li2(1/4)": Li2(mp.mpf(1)/4), "Li2(-1/3)": Li2(-mp.mpf(1)/3),
     "Li2(1/9)": Li2(mp.mpf(1)/9), "Li2(-1/8)": Li2(-mp.mpf(1)/8)}
B = {"Li2(2r2/3)": Li2(2*r2/3), "Li2(-b^4)": Li2(-b**4), "Li2(b^2/3)": Li2(b**2/3), "Li2(-b^2/2)x": None}
B = {k: v for k, v in B.items() if v is not None}
best = []
for an, av in A.items():
    for bn, bv in B.items():
        basis = prods + [("Li2(b)", Li2(b)), (an, av), (bn, bv)]
        names, vals = [], []
        for n, v in basis:
            names += [n, "r2*"+n]; vals += [v, r2*v]
        for tname, t in (("Q0", Q0), ("s2", s2)):
            out = find_relation(t, vals, names, 300, maxcoeff=10**7)
            if out:
                rel, res = out
                h = max(abs(c) for c in rel); nz = sum(1 for c in rel if c)
                best.append((tname, h, nz, an, bn, rel, names, res))
for tname in ("Q0", "s2"):
    cand = sorted([x for x in best if x[0] == tname], key=lambda x: (x[1], x[2]))
    for c in cand[:4]:
        tname, h, nz, an, bn, rel, names, res = c
        c0 = rel[0]
        terms = ["(%s)%s" % (mp.nstr(mp.mpf(-ci)/c0, 10) if (-ci) % c0 else str(-ci//c0), n) for ci, n in zip(rel[1:], names) if ci]
        print("%s  height=%d nonzero=%d basis[%s,%s] residual=%s\n   %s = (1/%d)*[%s]" % (tname, h, nz, an, bn, mp.nstr(res, 3), tname, c0,
              " + ".join("%d*%s" % (-ci, n) for ci, n in zip(rel[1:], names) if ci)))
