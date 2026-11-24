"""w6 identification step 19: c_3 against the large template family at D=850 (dps-900 run) and D=700 (same data truncated is NOT independent,
so the second precision uses the dps-400 run with a subfamily where possible).  Family groups (all with rational coefficients; sqrt2 variant optional):
  G   : Gamma-world L-family wt 6  (monomials in lam, L, pi^2, z3, z5)                           25
  L2  : Li4(1/2)*{lam^2, lam L, L^2, pi^2}, Li5(1/2)*{lam, L}                                       6
  A   : Li2(b) * wt-4 monomials in {lam, L, pi^2, z3}                                              11
  T3  : Li3(x) * wt-3 monomials, x in {b, -b, 1/sqrt2}                                              21
  T4  : Li4(x) * wt-2 monomials, x in {b,-b,1/sqrt2,2-sqrt2,b/sqrt2}                                20
  T5  : Li5(x) * {lam, L}, x in {b,-b,1/sqrt2,2-sqrt2,b/sqrt2}                                      10
  T6  : Li6(x), same x                                                                               5
  AA  : Li2(b)^2 * wt-2 monomials (4) ; Li2(b)*Li3(x)*{lam, L}, x in {b,-b,1/sqrt2} (6)           10
  P   : Li3(x)Li3(y) (6), Li2(b)*Li4(x) (5), Li2(b)^3 (1);   L3: Li6(1/2)
usage: s26_c3_bigfamily.py GROUPS RUN D     e.g.  G,T5 900 850
Printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel2
groups = sys.argv[1].split(","); RUN = int(sys.argv[2]); D = int(sys.argv[3])
mp.mp.dps = D + 30
print("mp.mp.dps =", mp.mp.dps, "groups", groups, flush=True)
def mono(w, gens):
    out = []
    def rec(i, wleft, name, val):
        if wleft == 0: out.append((name or "1", val)); return
        if i == len(gens): return
        nm, v, wt = gens[i]; e = 0
        while e * wt <= wleft:
            rec(i + 1, wleft - e * wt, name + (f"*{nm}^{e}" if e else ""), val * v ** e); e += 1
    rec(0, w, "", mp.mpf(1)); return out
r2 = mp.sqrt(2); b = r2 - 1; h = mp.mpf(1) / 2
lam = mp.log(2); L = mp.log(1 + r2); pi2 = mp.pi ** 2; z3 = mp.zeta(3); z5 = mp.zeta(5)
g4 = [("lam", lam, 1), ("L", L, 1), ("pi2", pi2, 2), ("z3", z3, 3), ("z5", z5, 5)]
A = mp.polylog(2, b)
args3 = {"b": b, "-b": -b, "1/r2": 1 / r2}
args5 = {"b": b, "-b": -b, "1/r2": 1 / r2, "2-r2": 2 - r2, "b/r2": b / r2}
B = []
if "G" in groups: B += mono(6, g4)
if "L2" in groups:
    Li4h = mp.polylog(4, h); Li5h = mp.polylog(5, h)
    B += [(f"Li4h*{n}", Li4h * v) for n, v in mono(2, g4[:3])] + [("Li5h*lam", Li5h * lam), ("Li5h*L", Li5h * L)]
if "A" in groups: B += [(f"A*{n}", A * v) for n, v in mono(4, g4[:4])]
if "T3" in groups:
    for nm, x in args3.items():
        li = mp.polylog(3, x); B += [(f"Li3({nm})*{n}", li * v) for n, v in mono(3, g4[:4])]
if "T4" in groups:
    for nm, x in args5.items():
        li = mp.polylog(4, x); B += [(f"Li4({nm})*{n}", li * v) for n, v in mono(2, g4[:3])]
if "T5" in groups:
    for nm, x in args5.items():
        li = mp.polylog(5, x); B += [(f"Li5({nm})*lam", li * lam), (f"Li5({nm})*L", li * L)]
if "T6" in groups:
    for nm, x in args5.items(): B.append((f"Li6({nm})", mp.polylog(6, x)))
if "AA" in groups:
    B += [(f"A^2*{n}", A ** 2 * v) for n, v in mono(2, g4[:3])]
    for nm, x in args3.items():
        li = mp.polylog(3, x); B += [(f"A*Li3({nm})*lam", A * li * lam), (f"A*Li3({nm})*L", A * li * L)]
if "P" in groups:
    l3 = {nm: mp.polylog(3, x) for nm, x in args3.items()}
    nms = list(l3)
    for ii in range(3):
        for jj in range(ii, 3):
            B.append((f"Li3({nms[ii]})*Li3({nms[jj]})", l3[nms[ii]] * l3[nms[jj]]))
    B += [(f"A*Li4({nm})", A * mp.polylog(4, x)) for nm, x in args5.items()]
    B += [("A^3", A ** 3)]
if "L3" in groups: B.append(("Li6h", mp.polylog(6, h)))
if "SQ" in groups:
    B += [("r2*" + n, r2 * v) for n, v in B]
cs = [mp.mpf(coef(RUN, j)) for j in (3, 2)]
mp.mp.dps = D + 30
n = len(B) + 1
print(f"family size {n} (incl. target), D={D}, run {RUN}, spurious size ~10^{D/n:.1f}", flush=True)
for j, cj in zip((3, 2), cs):
    good, internal = find_rel2(cj, B, D)
    print(f"c_{j}: {len(good)} relations with target, {len(internal)} internal relations")
    for r in good[:2]:
        print(f"   max|coeff|={r['maxcoeff']} residual {mp.nstr(r['res'],3)} target coeff {r['rel'][0]}")
        if r["maxcoeff"] < 10 ** 7: print("   TERMS:", r["terms"])
    for r in internal[:2]:
        print(f"   internal: max|coeff|={r['maxcoeff']} {r['terms'][:6]}...")
    if j == 3 and not good: pass
