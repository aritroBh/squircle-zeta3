"""w6 identification step 16 (template test at lower weight):
 (a) weight 3: c_0^{(3)} = W3' + 3S (closed form of thm:M3 in main.tex, transcribed) against the family
     { Gamma-world/L-family wt 3 , Li3(b), Li3(-b), Li3(1/sqrt2), Li2(b)*{lam, L} } x {1, sqrt2}  by LLL at D=380 and D=300.
 (b) weight 4: independence of the candidate family F4, then c_1^{(4)} and c_0^{(4)} (and the leading c_2^{(4)} as a control).
All printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel, lll_relations
r2f = lambda: mp.sqrt(2)
def setup(D):
    mp.mp.dps = D + 25
    S = {}
    r2 = mp.sqrt(2); b = r2 - 1; h = mp.mpf(1) / 2
    S.update(r2=r2, b=b, lam=mp.log(2), L=mp.log(1 + r2), pi2=mp.pi ** 2, z3=mp.zeta(3), z5=mp.zeta(5))
    S["A"] = mp.polylog(2, b)                               # Li2(b)
    S["B1"] = mp.polylog(3, b); S["B2"] = mp.polylog(3, -b); S["B3"] = mp.polylog(3, 1 / r2)
    S["Li4h"] = mp.polylog(4, h)
    return S
def mono(S, w, gens):
    """all monomials of weight w in the generators gens = list of (name, value, weight)"""
    out = []
    def rec(i, wleft, name, val):
        if wleft == 0:
            out.append((name or "1", val)); return
        if i == len(gens): return
        nm, v, wt = gens[i]
        e = 0
        while e * wt <= wleft:
            rec(i + 1, wleft - e * wt, name + (f"*{nm}^{e}" if e else ""), val * v ** e)
            e += 1
    rec(0, w, "", mp.mpf(1))
    return out
def with_sqrt2(B, S):
    return B + [("r2*" + n, S["r2"] * v) for n, v in B]
# ---------------- (a) weight 3
for D in (380, 300):
    S = setup(D)
    r2, b, lam, L, pi2, z3 = S["r2"], S["b"], S["lam"], S["L"], S["pi2"], S["z3"]
    S3 = z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12
    W3p = ((7 * r2 / 32 - mp.mpf(3) / 16) * z3 + (1 + r2) * S["B1"] + (1 - r2) * S["B2"] - 2 * r2 * S["B3"]
           + (2 * L - r2 * lam) * S["A"] + (r2 / 12 - mp.mpf(1) / 8) * pi2 * L - (1 + r2) / 24 * pi2 * lam
           + (1 + r2) / 24 * lam ** 3 + r2 / 4 * lam ** 2 * L - (1 + r2) / 2 * lam * L ** 2 + (5 - r2) / 6 * L ** 3)
    c0 = W3p + 3 * S3
    gens = [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3)]
    B = mono(S, 3, gens)
    B += [("A*lam", S["A"] * lam), ("A*L", S["A"] * L), ("Li3(b)", S["B1"]), ("Li3(-b)", S["B2"]), ("Li3(1/r2)", S["B3"])]
    B = with_sqrt2(B, S)
    print(f"[wt3] D={D}: family size {len(B)}, mp.mp.dps={mp.mp.dps}")
    R = find_rel(c0, B, D)
    if R:
        r = R[0]; print(f"   best relation: max|coeff|={r['maxcoeff']} residual {mp.nstr(r['res'],3)}  target coeff {r['rel'][0]}")
        if r["maxcoeff"] < 10 ** 5: print("   ", r["terms"])
    else:
        print("   no relation")
# ---------------- (b) weight 4
lc = json.load(open(os.path.join(W, "lower_coef_400.json")))["W4"]     # c_0, c_1, c_2
for D in (380, 300):
    S = setup(D)
    r2, b, lam, L, pi2, z3 = S["r2"], S["b"], S["lam"], S["L"], S["pi2"], S["z3"]
    gens = [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3), ("Li4h", S["Li4h"], 4)]
    B = mono(S, 4, gens)
    B += [(f"A*{n}", S["A"] * v) for n, v in mono(S, 2, gens[:3])]
    B += [("A^2", S["A"] ** 2)]
    for nm, key in (("Li3(b)", "B1"), ("Li3(-b)", "B2"), ("Li3(1/r2)", "B3")):
        B += [(f"{nm}*lam", S[key] * lam), (f"{nm}*L", S[key] * L)]
    x = {"b": b, "-b": -b, "1/r2": 1 / r2, "2-r2": 2 - r2, "b/r2": b / r2}
    for nm, xv in x.items():
        B.append((f"Li4({nm})", mp.polylog(4, xv)))
    n0 = len(B)
    B = with_sqrt2(B, S)
    print(f"[wt4] D={D}: family size {len(B)} (K-dim {n0}), spurious size ~10^{D/(len(B)+1):.1f}; mp.mp.dps={mp.mp.dps}")
    # independence of the family itself
    rows = lll_relations([v for _, v in B], D)
    nrel = sum(1 for rel, res, last in rows if abs(res) < mp.mpf(10) ** (-(D - 25)) and max(abs(a) for a in rel) < 10 ** 6)
    print(f"   relations inside the family itself (residual<1e-{D-25}, coeffs<1e6): {nrel}")
    names = [n for n, _ in B]
    for rel, res, last in rows:
        if abs(res) < mp.mpf(10) ** (-(D - 25)) and max(abs(a) for a in rel) < 10 ** 6:
            print("     internal:", [(names[i], a) for i, a in enumerate(rel) if a != 0])
    for j in (2, 1, 0):
        cj = mp.mpf(lc[j]); mp.mp.dps = D + 25
        R = find_rel(cj, B, D)
        if R:
            r = R[0]; print(f"   c_{j}^(4): best relation max|coeff|={r['maxcoeff']} residual {mp.nstr(r['res'],3)} target coeff {r['rel'][0]}")
            if r["maxcoeff"] < 10 ** 5: print("      ", r["terms"])
        else:
            print(f"   c_{j}^(4): no relation")
