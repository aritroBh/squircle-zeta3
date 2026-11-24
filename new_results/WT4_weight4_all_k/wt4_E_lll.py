"""WT4 step E: LLL closed forms.
D2: level-2 (alternating MZV) weight-4 basis {z3*log2, pi^4, pi^2 log2^2, log2^4, Li4(1/2)} (+sqrt2 twists as control).
D1: the full M4a 60-number basis (22 products x2 + 8 Li4 args x2).
Two digit-precisions each; accept only identical integer relations; flint lll gram='exact'."""
import os, sys, time, itertools, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
import flint
from fractions import Fraction as Fr

WORK = 930
mp.mp.dps = WORK
print("mp.mp.dps =", mp.mp.dps, flush=True)
t0 = time.time()

r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1
pi2 = mp.pi**2; Li2b = mp.polylog(2, b)
z3 = mp.zeta(3); Li3b = mp.polylog(3, b); Li3mb = mp.polylog(3, -b); Li3h = mp.polylog(3, 1/r2)
w1 = [("log2", lam), ("L", L)]
w2 = [("pi^2", pi2), ("Li2(b)", Li2b)]
w3 = [("z3", z3), ("Li3(b)", Li3b), ("Li3(-b)", Li3mb), ("Li3(1/sqrt2)", Li3h)]
prods = []
for c in itertools.combinations_with_replacement(range(2), 4):
    prods.append(("*".join(w1[i][0] for i in c), mp.fprod(w1[i][1] for i in c)))
for (n2, v2) in w2:
    for c in itertools.combinations_with_replacement(range(2), 2):
        prods.append((n2+"*"+"*".join(w1[i][0] for i in c), v2*mp.fprod(w1[i][1] for i in c)))
for c in itertools.combinations_with_replacement(range(2), 2):
    prods.append(("*".join(w2[i][0] for i in c), mp.fprod(w2[i][1] for i in c)))
for (n3, v3) in w3:
    for (n1, v1) in w1:
        prods.append((n3+"*"+n1, v3*v1))
print(f"{len(prods)} products; t={time.time()-t0:.1f}", flush=True)

def q2(p, q="0"):
    return mp.mpf(p) + r2*mp.mpf(q)
KEPT = [("Li4(1/2)", q2("1/2")), ("Li4(b)", q2("-1","1")), ("Li4(-b)", q2("1","-1")),
        ("Li4(1/sqrt2)", q2("0","1/2")), ("Li4(-b^2)", q2("-3","2")), ("Li4(2-sqrt2)", q2("2","-1")),
        ("Li4(1-1/sqrt2)", q2("1","-1/2")), ("Li4(-1+1/sqrt2)", q2("-1","1/2"))]
kept = [(nm, mp.polylog(4, z)) for nm, z in KEPT]
print(f"8 Li4 values; t={time.time()-t0:.1f}", flush=True)
Li4half = kept[0][1]

def qs2(lst):
    out = []
    for nm, v in lst:
        out += [(nm, v), ("sqrt2*"+nm, r2*v)]
    return out

def lll(vals, digits):
    C = mp.mpf(10)**digits
    n = len(vals)
    rows = [[1 if i == j else 0 for j in range(n)] + [int(mp.nint(C*vals[i]))] for i in range(n)]
    Bm = flint.fmpz_mat(rows).lll(gram='exact')
    best = [int(Bm[0, j]) for j in range(n)]
    with mp.workdps(digits + 30):
        resid = mp.fsum(c*v for c, v in zip(best, vals))
    return best, max(abs(c) for c in best), resid

D1V = mp.mpf(open(os.path.join(HERE, "d1_860.txt")).read().strip())
D2V = mp.mpf(open(os.path.join(HERE, "d2_860.txt")).read().strip())

out = {}

# ---------- D2: level-2 basis (5 numbers) ----------
lvl2 = [("z3*log2", z3*lam), ("pi^4", pi2**2), ("pi^2*log2^2", pi2*lam**2),
        ("log2^4", lam**4), ("Li4(1/2)", Li4half)]
for D in (700, 780):
    full = [("D2", D2V)] + lvl2
    vec, mx, res = lll([v for _, v in full], D)
    print(f"[D2 lvl2] LLL at {D} digits: D2 coef={vec[0]}, max|coef|={mx:.3e}, residual={mp.nstr(res,3)} (mp.mp.dps={mp.mp.dps})", flush=True)
    out.setdefault("D2_lvl2", []).append({"D": D, "vec": vec, "maxcoef": str(mx), "residual": mp.nstr(res, 6)})
v1, v2 = out["D2_lvl2"][0]["vec"], out["D2_lvl2"][1]["vec"]
same = (v1 == v2) or (v1 == [-x for x in v2])
print(f"[D2 lvl2] identical at 700/780: {same}", flush=True)
assert same and v1[0] != 0
g = 0
for x in v1: g = __import__('math').gcd(g, abs(x))
rel = {nm: Fr(-c, v1[0]).__reduce__() and Fr(-c, v1[0]) for c, nm in zip(v1[1:], [n for n,_ in lvl2]) if c}
# normalize: D2 = sum Fr * basis
rel = {nm: Fr(-c, v1[0]) for c, nm in zip(v1[1:], [n for n,_ in lvl2]) if c}
out["D2_closed_form"] = {nm: str(f) for nm, f in rel.items()}
with mp.workdps(880):
    chk = D2V - mp.fsum(f*v for (nm, v), f in zip(lvl2, [Fr(out["D2_closed_form"].get(nm, "0")) for nm,_ in lvl2]))
print(f"[D2 lvl2] closed form - D2(860d) = {mp.nstr(chk,3)} at workdps 880", flush=True)
out["D2_check_860"] = mp.nstr(chk, 6)

# ---------- D2 control: same basis x {1,sqrt2} (10 numbers) ----------
for D in (700, 780):
    full = [("D2", D2V)] + qs2(lvl2)
    vec, mx, res = lll([v for _, v in full], D)
    print(f"[D2 lvl2x2] LLL at {D}: D2 coef={vec[0]}, max|coef|={mx:.3e}, residual={mp.nstr(res,3)}", flush=True)
    out.setdefault("D2_lvl2x2", []).append({"D": D, "vec": vec, "maxcoef": str(mx), "residual": mp.nstr(res, 6)})

# ---------- D1: full 60-number basis ----------
full_basis = qs2(prods) + qs2(kept)
print(f"D1 basis size = {1+len(full_basis)} (expect 61); t={time.time()-t0:.1f}", flush=True)
for D in (700, 780):
    full = [("D1", D1V)] + full_basis
    t1 = time.time()
    vec, mx, res = lll([v for _, v in full], D)
    print(f"[D1 full60] LLL at {D} digits ({time.time()-t1:.0f}s): D1 coef={vec[0]}, max|coef|={mx:.3e}, residual={mp.nstr(res,3)}", flush=True)
    print("   relation:", " ".join(f"{c:+d}*{nm}" for c, nm in zip(vec, [n for n,_ in full]) if c), flush=True)
    out.setdefault("D1_full60", []).append({"D": D, "vec": vec, "maxcoef": str(mx), "residual": mp.nstr(res, 6)})
v1, v2 = out["D1_full60"][0]["vec"], out["D1_full60"][1]["vec"]
same = (v1 == v2) or (v1 == [-x for x in v2])
print(f"[D1 full60] identical at 700/780: {same}", flush=True)
assert same and v1[0] != 0
names = [n for n,_ in full_basis]
rel1 = {nm: Fr(-c, v1[0]) for c, nm in zip(v1[1:], names) if c}
out["D1_closed_form"] = {nm: str(f) for nm, f in rel1.items()}
with mp.workdps(880):
    chk = D1V - mp.fsum(Fr(out["D1_closed_form"].get(nm, "0"))*v for nm, v in full_basis)
print(f"[D1 full60] closed form - D1(860d) = {mp.nstr(chk,3)} at workdps 880", flush=True)
out["D1_check_860"] = mp.nstr(chk, 6)

json.dump(out, open(os.path.join(HERE, "lll_results.json"), "w"), indent=1)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_E lll: D2 lvl2 identical={same}; D2 chk={out['D2_check_860']}; D1 chk={out['D1_check_860']}\n")
print("DONE", time.time()-t0, "s")
