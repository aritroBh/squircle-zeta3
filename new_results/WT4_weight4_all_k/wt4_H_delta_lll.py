"""WT4 step H (bonus): LLL identification of the k=3 offset delta = P(3) - W_4(3)
in the M4a 60-number weight-4 basis, at two digit-precisions (identical relation required)."""
import os, sys, time, itertools, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
import flint
from fractions import Fraction as Fr

mp.mp.dps = 920
print("mp.mp.dps =", mp.mp.dps, flush=True)
t0 = time.time()
r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1
pi2 = mp.pi**2; z3 = mp.zeta(3)
Li2b = mp.polylog(2, b); Li3b = mp.polylog(3, b); Li3mb = mp.polylog(3, -b); Li3h = mp.polylog(3, 1/r2)
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
def q2(p, q="0"):
    return mp.mpf(p) + r2*mp.mpf(q)
KEPT = [("Li4(1/2)", q2("1/2")), ("Li4(b)", q2("-1","1")), ("Li4(-b)", q2("1","-1")),
        ("Li4(1/sqrt2)", q2("0","1/2")), ("Li4(-b^2)", q2("-3","2")), ("Li4(2-sqrt2)", q2("2","-1")),
        ("Li4(1-1/sqrt2)", q2("1","-1/2")), ("Li4(-1+1/sqrt2)", q2("-1","1/2"))]
kept = [(nm, mp.polylog(4, z)) for nm, z in KEPT]
def qs2(lst):
    out = []
    for nm, v in lst:
        out += [(nm, v), ("sqrt2*"+nm, r2*v)]
    return out
full_basis = qs2(prods) + qs2(kept)

# delta at ~810 digits: P(3) from files; W_4(3) from thm:W4 closed form at 920 dps
W4pp = mp.mpf(open("<project>/research/08_theory/M4a_decomposition/W4pp_1060.txt").read().strip())
D1V = mp.mpf(open(os.path.join(HERE, "d1_860.txt")).read().strip())
D2V = mp.mpf(open(os.path.join(HERE, "d2_860.txt")).read().strip())
import w4_symbol as ws
def atom_eval(a):
    if a[0] == 'li':
        z = a[2]; return mp.polylog(a[1], mp.mpf(str(z.p)) + r2*mp.mpf(str(z.q)))
    if a[0] == 'log':
        j2, m = ws.logvec(a[1]); return mp.mpf(float(j2))*lam + mp.mpf(float(m))*L
    if a == ('pi2',): return pi2
    if a == ('z3',): return z3
W43 = mp.mpf(0)
for c, atoms in ws.conjecture_terms():
    cm = mp.mpf(str(c.p)) + r2*mp.mpf(str(c.q))
    v = cm
    for a in atoms: v *= atom_eval(a)
    W43 += v
delta = (W4pp - D1V + D2V) - W43
print(f"delta = {mp.nstr(delta, 60)}  (t={time.time()-t0:.0f}s)", flush=True)

def lll(vals, digits):
    C = mp.mpf(10)**digits
    n = len(vals)
    rows = [[1 if i == j else 0 for j in range(n)] + [int(mp.nint(C*vals[i]))] for i in range(n)]
    Bm = flint.fmpz_mat(rows).lll(gram='exact')
    best = [int(Bm[0, j]) for j in range(n)]
    with mp.workdps(digits + 30):
        resid = mp.fsum(c*v for c, v in zip(best, vals))
    return best, max(abs(c) for c in best), resid

out = {"delta_60d": mp.nstr(delta, 60)}
vecs = []
for D in (700, 780):
    t1 = time.time()
    full = [("delta", delta)] + full_basis
    vec, mx, res = lll([v for _, v in full], D)
    vecs.append(vec)
    print(f"[delta full60] LLL at {D} ({time.time()-t1:.0f}s): delta coef={vec[0]}, max|coef|={mx:.3e}, residual={mp.nstr(res,3)}", flush=True)
    out.setdefault("runs", []).append({"D": D, "vec": vec, "maxcoef": str(mx), "residual": mp.nstr(res, 6)})
same = (vecs[0] == vecs[1]) or (vecs[0] == [-x for x in vecs[1]])
print("identical at 700/780:", same, flush=True)
if same and vecs[0][0] != 0:
    names = [n for n, _ in full_basis]
    rel = {nm: Fr(-c, vecs[0][0]) for c, nm in zip(vecs[0][1:], names) if c}
    out["delta_closed_form_candidate"] = {nm: str(f) for nm, f in rel.items()}
    print(f"candidate: {len(rel)} terms; max |coef| = {max(abs(f.numerator) for f in rel.values())}")
    with mp.workdps(850):
        chk = delta - mp.fsum(f*v for nm, v in ((nm, v) for nm, v in full_basis) if (f := rel.get(nm)))
    # safer recompute without walrus-in-genexp over dict
    with mp.workdps(850):
        s = mp.mpf(0)
        for nm, v in full_basis:
            if nm in rel: s += rel[nm]*v
        chk = delta - s
    print(f"candidate - delta = {mp.nstr(chk,3)} at workdps 850", flush=True)
    out["check_810"] = mp.nstr(chk, 6)
json.dump(out, open(os.path.join(HERE, "delta_k3_lll.json"), "w"), indent=1)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_H delta LLL: identical={same}, terms={len(out.get('delta_closed_form_candidate', {}))}, check={out.get('check_810')}\n")
print("DONE", time.time()-t0, "s")
