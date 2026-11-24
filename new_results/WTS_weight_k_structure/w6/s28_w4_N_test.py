"""w6 identification step 21: which of the 8 weight-4 basic integrals N^(4)_{key,r} lie in the (K-closed) 56-element family used for c_1^(4)?
D=380.  Then LLL on all 8 N together with the family to count relations (dimension of the span modulo the family)."""
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
for nm, xv in {"b": b, "-b": -b, "1/r2": 1 / r2, "2-r2": 2 - r2, "b/r2": b / r2}.items():
    B.append((f"Li4({nm})", mp.polylog(4, xv)))
B = with_sqrt2(B, S)
N = json.load(open(os.path.join(W, "Nints4_400.json")))["N"]
mp.mp.dps = D + 25
print("family size", len(B), " mp.mp.dps", mp.mp.dps)
for k, v in N.items():
    good, internal = find_rel2(mp.mpf(v), B, D)
    print(f"N4[{k}]: " + (f"relation max|coeff|={good[0]['maxcoeff']} res {mp.nstr(good[0]['res'],3)}" if good and good[0]['maxcoeff'] < 10**5 else "NOT in family"))
# joint: all 8 N + family
vals = [mp.mpf(v) for v in N.values()] + [v for _, v in B]
names = [f"N4[{k}]" for k in N] + [n for n, _ in B]
rows = lll_relations(vals, D)
print("joint relations (coeff < 1e5):")
for rel, res, last in rows:
    if abs(res) < mp.mpf(10) ** (-(D - 25)) and max(abs(a) for a in rel) < 10 ** 5:
        print("  ", [(names[i], a) for i, a in enumerate(rel) if a != 0][:30])
