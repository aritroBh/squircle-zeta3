"""w6 identification step 23: weight 4 with the extra Li4 argument -b^2 (the 7th depth-1 weight-4 generator: l_4 = 7 = {1/2, b, -b, 1/sqrt2, 2-sqrt2, b/sqrt2, -b^2}).
Counts the relations among the 8 weight-4 basic integrals N and the K-closed family (57 numbers) and tests c_0^(4), c_1^(4).  D=380."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel2, lll_relations
exec(open(os.path.join(W, "s23_template.py")).read().split("# ---------------- (a) weight 3")[0].split("from w6lll import find_rel, lll_relations")[1])
D = 380
S = setup(D)
r2, b, lam, L, pi2, z3 = S["r2"], S["b"], S["lam"], S["L"], S["pi2"], S["z3"]
gens = [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3), ("Li4h", S["Li4h"], 4)]
B = mono(S, 4, gens)
B += [(f"A*{n}", S["A"] * v) for n, v in mono(S, 2, gens[:3])] + [("A^2", S["A"] ** 2)]
for nm, key in (("Li3(b)", "B1"), ("Li3(-b)", "B2"), ("Li3(1/r2)", "B3")):
    B += [(f"{nm}*lam", S[key] * lam), (f"{nm}*L", S[key] * L)]
for nm, xv in {"b": b, "-b": -b, "1/r2": 1 / r2, "2-r2": 2 - r2, "b/r2": b / r2, "-b^2": -b * b}.items():
    B.append((f"Li4({nm})", mp.polylog(4, xv)))
B = with_sqrt2(B, S)
N = json.load(open(os.path.join(W, "Nints4_400.json")))["N"]
lc = json.load(open(os.path.join(W, "lower_coef_400.json")))["W4"]
mp.mp.dps = D + 25
print("family size", len(B), "spurious size ~10^%.1f" % (D / (len(B) + 8)))
vals = [mp.mpf(v) for v in N.values()] + [v for _, v in B]
names = [f"N4[{k}]" for k in N] + [n for n, _ in B]
rows = lll_relations(vals, D)
cred = [(rel, res) for rel, res, last in rows if abs(res) < mp.mpf(10) ** (-(D - 25)) and max(abs(a) for a in rel) < 10 ** 4.5]
print("credible joint relations:", len(cred))
for rel, res in cred:
    print("  max", max(abs(a) for a in rel), [(names[i], a) for i, a in enumerate(rel) if a != 0 and i < 8], "| family terms:", sum(1 for i, a in enumerate(rel) if a != 0 and i >= 8))
for j in (1, 0):
    cj = mp.mpf(lc[j]); mp.mp.dps = D + 25
    good, internal = find_rel2(cj, B, D)
    print(f"c_{j}^(4): " + (f"relation max|coeff|={good[0]['maxcoeff']} res {mp.nstr(good[0]['res'],3)} (internal relations: {len(internal)})" if good else "none"))
    if good and good[0]["maxcoeff"] < 10 ** 5:
        print("   TERMS:", good[0]["terms"])
