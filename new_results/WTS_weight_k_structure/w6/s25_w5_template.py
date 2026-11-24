"""w6 identification step 18: weight-5 next-to-leading coefficient c_2^(5) (and leading c_3^(5) as control) against growing families, to learn
the structure of the subleading coefficient before attacking weight 6.  LLL at D=380 (dps-400 data) and D=300.  Printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel
lc = json.load(open(os.path.join(W, "lower_coef_400.json")))["W5"]      # c_0..c_3 of W_5(k), k>=5
def mono(w, gens):
    out = []
    def rec(i, wleft, name, val):
        if wleft == 0: out.append((name or "1", val)); return
        if i == len(gens): return
        nm, v, wt = gens[i]; e = 0
        while e * wt <= wleft:
            rec(i + 1, wleft - e * wt, name + (f"*{nm}^{e}" if e else ""), val * v ** e); e += 1
    rec(0, w, "", mp.mpf(1)); return out
def consts_w5(D):
    mp.mp.dps = D + 25
    r2 = mp.sqrt(2); b = r2 - 1; h = mp.mpf(1) / 2
    C = dict(r2=r2, b=b, lam=mp.log(2), L=mp.log(1 + r2), pi2=mp.pi ** 2, z3=mp.zeta(3), z5=mp.zeta(5),
             Li4h=mp.polylog(4, h), Li5h=mp.polylog(5, h), A=mp.polylog(2, b))
    C["args"] = {"b": b, "-b": -b, "1/r2": 1 / r2, "2-r2": 2 - r2, "b/r2": b / r2}
    return C
def fam(D, groups):
    C = consts_w5(D)
    lam, L, pi2, z3, z5 = C["lam"], C["L"], C["pi2"], C["z3"], C["z5"]
    B = []
    if "gamma" in groups:                     # Gamma-world L-family, weight 5 (monomials in lam, L, pi2, z3, z5)
        B += mono(5, [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3), ("z5", z5, 5)])
    if "lvl2" in groups:                      # non-Gamma level-2 weight-5 pieces
        B += [("Li4h*lam", C["Li4h"] * lam), ("Li4h*L", C["Li4h"] * L), ("Li5h", C["Li5h"])]
    if "li4" in groups:
        for nm, x in C["args"].items():
            v = mp.polylog(4, x); B += [(f"Li4({nm})*lam", v * lam), (f"Li4({nm})*L", v * L)]
    if "li3" in groups:
        for nm, x in list(C["args"].items())[:3]:
            v = mp.polylog(3, x)
            for n2, v2 in mono(2, [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2)]): B.append((f"Li3({nm})*{n2}", v * v2))
    if "li5" in groups:
        for nm, x in list(C["args"].items())[:3]: B.append((f"Li5({nm})", mp.polylog(5, x)))
    if "A" in groups:
        for n2, v2 in mono(3, [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3)]): B.append((f"A*{n2}", C["A"] * v2))
    if "sq" in groups:
        B += [("r2*" + n, C["r2"] * v) for n, v in B]
    return B
for groups in (["gamma", "li4"], ["gamma", "li4", "li3"], ["gamma", "li4", "li3", "lvl2"], ["gamma", "li4", "li3", "lvl2", "li5"],
               ["gamma", "li4", "li3", "lvl2", "li5", "A"]):
    for j in (3, 2):
        for D in (380, 300):
            B = fam(D, groups); cj = mp.mpf(lc[j]); mp.mp.dps = D + 25
            R = find_rel(cj, B, D); n = len(B) + 1
            msg = f"groups={groups} c_{j}^(5) D={D} n={n} spurious~10^{D/n:.1f}: "
            msg += "no relation" if not R else f"best max|coeff|={R[0]['maxcoeff']} res {mp.nstr(R[0]['res'],3)} target coeff {R[0]['rel'][0]}"
            print(msg, flush=True)
            if R and R[0]["maxcoeff"] < 10 ** 5 and D == 380: print("    TERMS:", R[0]["terms"], flush=True)
