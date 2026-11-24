# 22_lll_vertex_try.py -- exploratory LLL search for the unidentified vertex integrals V1, V3
# (and the combination U = sqrt2 (V3 - V2) + 2 V1 entering Q1) in a guessed weight<=3 basis.
import mpmath as mp, json, sys, itertools
from lll_util import lll_relation, describe
fn = sys.argv[1] if len(sys.argv) > 1 else "Q1_sub_400.json"
if fn.endswith("1250"):
    d = {"V1": json.load(open("V1_1250.json"))["V1"], "V3p": json.load(open("V3p_1250.json"))["V3p"]}
else:
    d0 = json.load(open(fn)); d = {"V1": d0["V1"], "V3p": str(mp.mpf(d0["V3"]) - mp.mpf(d0["V2"])) if False else None}
    mp.mp.dps = 450; d["V3p"] = mp.nstr(mp.mpf(d0["V3"]) - mp.mpf(d0["V2"]), 400)
DPS = len(d["V1"]) - 5
mp.mp.dps = DPS + 20
r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1; al = mp.atan(1/r2)
Li2 = lambda x: mp.polylog(2, x); Li3 = lambda x: mp.polylog(3, x)
w = 1 - 1j*r2
Ti2 = mp.im(Li2(1j*r2))
RLi3w = mp.re(Li3(w)); ILi3w = mp.im(Li3(w))
w1 = [("lam",lam),("l3",l3),("L",L),("pi",pi),("al",al)]
w2 = [(x[0]+"*"+y[0], x[1]*y[1]) for x, y in itertools.combinations_with_replacement(w1, 2)] + \
     [("Li2b",Li2(b)),("Li2mb4",Li2(-b**4)),("Li2_13",Li2(mp.mpf(1)/3)),("Ti2r2",Ti2)]
w3 = [(x[0]+"*"+y[0]+"*"+z[0], x[1]*y[1]*z[1]) for x, y, z in itertools.combinations_with_replacement(w1, 3)] + \
     [(x[0]+"*"+y[0], x[1]*y[1]) for x in w1 for y in w2[-4:]] + \
     [("z3",mp.zeta(3)),("Li3b",Li3(b)),("Li3mb",Li3(-b)),("Li3r",Li3(1/r2)),("Li3mb4",Li3(-b**4)),("Li3_13",Li3(mp.mpf(1)/3)),
      ("ReLi3w",RLi3w),("ImLi3w",ILi3w)]
base = [("1",mp.mpf(1))] + w1 + w2 + w3
mode = sys.argv[2] if len(sys.argv) > 2 else "both"
names, vals = [], []
for n, v in base:
    if mode in ("both", "one"): names.append(n); vals.append(v)
    if mode in ("both", "r2"): names.append("r2*"+n); vals.append(r2*v)
V1 = mp.mpf(d["V1"]); V3p = mp.mpf(d["V3p"])
U = r2*V3p + 2*V1
print("basis size", len(vals), " digits", DPS)
for key, t in (("V1", V1), ("V3p", V3p), ("U", U)):
    rel, res, h = lll_relation([t] + vals, DPS - 10)
    expect = (DPS - 10)/(len(vals) + 1)
    sig = mp.log10(h) < 0.6*expect
    print("%s: height 10^%.2f  (random expectation ~10^%.1f)  -> %s" % (key, float(mp.log10(h)), float(expect),
          "SIGNIFICANT" if sig else "no relation (not significant)"))
    if sig:
        print("   ", describe(rel, names))
