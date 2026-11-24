import os
work='<project>/research/08_theory/WT4_weight4_all_k'
script = r'''"""WT4 standalone verifier. Reruns the controls and prints PASS/FAIL per check:
 [1] k=4 degree-3 kernels == M4a omega4.json (exact); quadratic kernel law spot-check k=7..10 (exact)
 [2] D1/D2/W4pp G-expressions rebuilt from deltas/kernels == d1.tsv/d2.tsv/gpl_w4pp.tsv (exact)
 [3] TSV vs mpmath quadrature at dps 60 (flint G-evaluator ctx.prec=500)
 [4] proved closed forms R_D1, R_D2 (200 dps) vs quadrature values (200 dps)
 [5] exact certificates: subprocess verify_wt4_cert.py D2 / D1 -> VERDICT PROVED
 [6] k=3 offset delta = W4''-D1+D2-W_4(3): nonzero (~0.15446753...), W_4(3) quad == thm:W4 eval
"""
import sys, os, json, subprocess, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
PY = sys.executable
RESULTS = []
def report(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)

t0 = time.time()
# ---------------- [1] exact kernel controls ----------------
import pickle
from fractions import Fraction as Fr
import fastred
from fastred import QSI, ZERO, AQ, BQ, reduce_k_fast, pf_zero, pf_add, pf_scale, pf_sub
POLE_NAME = {ZERO: "0", AQ: "a", -AQ: "-a", BQ: "b", -BQ: "-b"}
def pf_to_omega(P):
    real, ipos, ineg, const = P
    assert not ipos and not ineg and const.is_zero()
    out = []
    for r, c in real.items():
        assert c.is_real()
        out.append((POLE_NAME[r], str(c.p0), str(c.p1)))
    return sorted(out)
def reduce_deg3(k):
    ipf = pickle.load(open(os.path.join(HERE, f"binteg_{k}.pkl"), "rb"))
    res = reduce_k_fast(k, integrand_pf=ipf)
    return {key: res[key] for key in [(2, 1), (3, 0)]}
M4A = "<project>/research/08_theory/M4a_decomposition"
om4 = json.load(open(os.path.join(M4A, "omega4.json")))
r4 = reduce_deg3(4)
ok = all(pf_to_omega(r4[key]) == sorted((t[0], t[1][0], t[1][1]) for t in om4["%d,%d" % key]) for key in [(2,1),(3,0)])
report("control k=4 kernels == M4a omega4.json (exact)", ok)
deltas = json.load(open(os.path.join(HERE, "deltas_deg3.json")))
def K_om(kk, key):
    return sorted((t[0], t[1][0], t[1][1]) for t in json.load(open(os.path.join(HERE,"kernels_deg3.json")))[str(kk)][key])
def pf_from_om(terms):
    NAME = {"0": ZERO, "a": AQ, "-a": -AQ, "b": BQ, "-b": -BQ}
    P = pf_zero()
    for nm, p, q in terms:
        P[0][NAME[nm]] = P[0].get(NAME[nm], ZERO) + QSI(Fr(p), Fr(q))
    return P
D1pf = {key: pf_from_om([(nm, p, q) for nm, (p, q) in [(x[0], x[1]) for x in deltas["%d,%d"%key]["D1"]]]) for key in [(2,1),(3,0)]}
D2pf = {key: pf_from_om([(nm, p, q) for nm, (p, q) in [(x[0], x[1]) for x in deltas["%d,%d"%key]["D2"]]]) for key in [(2,1),(3,0)]}
K4pf = {key: pf_from_om(K_om(4, key)) for key in [(2,1),(3,0)]}
from math import comb
lawok = True
for k in range(7, 11):
    rk = reduce_deg3(k)
    for key in [(2,1),(3,0)]:
        pred = pf_add(pf_add(K4pf[key], pf_scale(D1pf[key], Fr(k-4))), pf_scale(D2pf[key], Fr(comb(k-4,2))))
        diff = pf_sub(pred, rk[key])
        lawok = lawok and (not diff[0]) and (not diff[1]) and (not diff[2]) and diff[3].is_zero()
report("quadratic kernel law k=7..10 (exact; k=11..14 done in wt4_A_kernels.py)", lawok)

# ---------------- [2] G-expression rebuild == TSVs ----------------
from collections import defaultdict
from qsqrt2 import QS2 as Q, A, B
from w4_symbol import log_product, LT, MT, clean
POLE = {"0": Q(0), "a": A, "-a": -A, "b": B, "-b": -B}
def g_expression(i, j, terms):
    out = defaultdict(lambda: Q(0))
    for name, (p, q) in terms:
        z = POLE[name]; coefficient = Q(Fr(p), Fr(q))
        for (w, k), d in log_product(*((LT,)*i + (MT,)*j)).items():
            out[((z,) + w, k)] += coefficient * d
    return clean(out)
def build(kerns):
    out = defaultdict(lambda: Q(0))
    for (i, j) in ((2, 1), (3, 0)):
        for kk, v in g_expression(i, j, kerns["%d,%d" % (i, j)]).items():
            out[kk] += v
    return clean(out)
import csv
def load_tsv_dict(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    out = defaultdict(lambda: Q(0))
    for r in rows:
        # word string equality via re-format: rebuild QS2 tuple string
        out[(r["word"], int(r["k"]))] += Q(Fr(r["coefficient_p"]), Fr(r["coefficient_q"]))
    return dict(out)
kernels = json.load(open(os.path.join(HERE, "kernels_deg3.json")))
gD1 = build({k: deltas[k]["D1"] for k in ["2,1","3,0"]})
gD2 = build({k: deltas[k]["D2"] for k in ["2,1","3,0"]})
gW4 = build(kernels["4"])
ok1 = {(str(w), k): c for (w, k), c in gD1.items()} == load_tsv_dict(os.path.join(HERE, "d1.tsv"))
ok2 = {(str(w), k): c for (w, k), c in gD2.items()} == load_tsv_dict(os.path.join(HERE, "d2.tsv"))
ok3 = {(str(w), k): c for (w, k), c in gW4.items()} == load_tsv_dict(os.path.join(M4A, "gpl_w4pp.tsv"))
report("G-expr(D1) == d1.tsv (exact)", ok1)
report("G-expr(D2) == d2.tsv (exact)", ok2)
report("G-expr(W4pp) == gpl_w4pp.tsv (exact control)", ok3)

# ---------------- [3] TSV vs quadrature at dps 60 ----------------
from flint import acb, acb_series, acb_poly, ctx
import mpmath as mp
from functools import cache
from fractions import Fraction as F
ctx.prec = 500; ctx.cap = 300; N = ctx.cap
mp.mp.dps = 60
def qa(c): return acb(str(c.p)) + acb(2).sqrt()*acb(str(c.q))
@cache
def shd(u, v):
    if not u: return ((v, 1),)
    if not v: return ((u, 1),)
    d = {}
    for w, m in shd(u[1:], v): d[(u[0],)+w] = d.get((u[0],)+w, 0) + m
    for w, m in shd(u, v[1:]): d[(v[0],)+w] = d.get((v[0],)+w, 0) + m
    return tuple(d.items())
@cache
def reg_words(w, e):
    if not w or w[0] != e: return ((w, F(1)),)
    terms = dict(shd((e,), w[1:])); cw = terms.pop(w); d = {}
    for u, m in terms.items():
        for t, cf in reg_words(u, e): d[t] = d.get(t, F(0)) + F(-m, cw)*cf
    return tuple((t, cf) for t, cf in d.items() if cf)
@cache
def real_series(poles):
    if not poles: return acb_series([1])
    tail = real_series(poles[1:])
    if poles[0] == Q(0): return acb_series([0]+[tail[j]/j for j in range(1, N)])
    return (tail/acb_series([-qa(poles[0]), 1])).integral()
@cache
def real_g_raw(w, c, e):
    if not w: return acb(1)
    if c == e: return acb(0)
    def taylor(w, c, e):
        if w[-1] == c: return None
        nz = [abs((p-c)/(e-c)) for p in w if p != c]
        if nz and min(nz) <= Q(F(3, 2)): return None
        poles = tuple((p-c)/(e-c) for p in w)
        return acb_poly(real_series(poles).coeffs())(acb(1))
    v = taylor(w, c, e)
    if v is not None: return v
    v = taylor(tuple(reversed(w)), e, c)
    if v is not None: return (-1)**len(w)*v
    mid = (c+e)/2
    return sum((real_g_raw(w[:k], mid, e)*real_g_raw(w[k:], c, mid) for k in range(len(w)+1)), acb(0))
def real_g(w, c, e):
    return sum((acb(cf.numerator)/cf.denominator*real_g_raw(u, c, e) for u, cf in reg_words(w, e)), acb(0))
def parse_word(s):
    vals = {"0": Q(0), "(1+1*sqrt2)": A, "(-1-1*sqrt2)": -A, "(-1+1*sqrt2)": B, "(1-1*sqrt2)": -B}
    s = s.strip(); body = s[1:-1].strip()
    if not body: return ()
    parts, depth, cur = [], 0, ""
    for ch in body:
        if ch == "(": depth += 1
        if ch == ")": depth -= 1
        if ch == "," and depth == 0: parts.append(cur.strip()); cur = ""
        else: cur += ch
    if cur.strip(): parts.append(cur.strip())
    return tuple(vals[p] for p in parts)
def evaluate(path):
    ONE = Q(1); log2 = acb(2).log(); tot = acb(0)
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    for r in rows:
        w = parse_word(r["word"]); k = int(r["k"]); c = Q(F(r["coefficient_p"]), F(r["coefficient_q"]))
        tot += qa(c)*(log2**k)*(real_g(w, ONE, A) - real_g(w, ONE, B))
    return tot
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
def quad_value(kerns):
    tot = mp.mpf(0)
    for key, i, j in [("2,1", 2, 1), ("3,0", 3, 0)]:
        ts = [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in kerns[key]]
        def f(v, ts=ts, i=i, j=j):
            l1 = mp.log((v - bN)*(v + aN)/(4*v))
            l2 = mp.log((aN - v)*(v + bN)/(4*v))
            return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
        tot += mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], maxdegree=12)
    return tot
for name, tsv, kerns in [("D1", "d1.tsv", {k: deltas[k]["D1"] for k in ["2,1","3,0"]}),
                         ("D2", "d2.tsv", {k: deltas[k]["D2"] for k in ["2,1","3,0"]}),
                         ("W4pp(ctl)", os.path.join(M4A, "gpl_w4pp.tsv"), kernels["4"])]:
    gv = evaluate(os.path.join(HERE, tsv) if not tsv.startswith("/") else tsv)
    qv = quad_value(kerns)
    got = mp.mpf(gv.real.mid().str(70, radius=False, more=True))
    diff = abs(got - qv)
    report(f"TSV-vs-quadrature {name} (dps {mp.mp.dps})", diff < mp.mpf('1e-45'), f"|G-quad|={mp.nstr(diff,3)} flint_rad={gv.real.rad().str(4)}")

# ---------------- [4] closed forms vs quadrature at 200 dps ----------------
with mp.workdps(200):
    rr2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+rr2); b = rr2-1
    pi2 = mp.pi**2; z3 = mp.zeta(3)
    bv = {"log2": lam, "L": L, "pi^2": pi2, "Li2(b)": mp.polylog(2, b), "z3": z3,
          "Li3(b)": mp.polylog(3, b), "Li3(-b)": mp.polylog(3, -b), "Li3(1/sqrt2)": mp.polylog(3, 1/rr2)}
    def ev(nm):
        if nm == "pi^4": return pi2**2
        if nm == "pi^2*log2^2": return pi2*lam**2
        if nm == "log2^4": return lam**4
        v = mp.mpf(1)
        for p in nm.split("*"): v *= bv[p]
        return v
    lll = json.load(open(os.path.join(HERE, "lll_results.json")))
    qD1 = quad_value({k: deltas[k]["D1"] for k in ["2,1","3,0"]})
    qD2 = quad_value({k: deltas[k]["D2"] for k in ["2,1","3,0"]})
    for nm_, cf_, qv_ in [("D1", lll["D1_closed_form"], qD1), ("D2", lll["D2_closed_form"], qD2)]:
        rv = mp.mpf(0)
        for nm, fs in cf_.items():
            f = Fr(fs)
            rv += f*(rr2*ev(nm[6:]) if nm.startswith("sqrt2*") else ev(nm))
        e = abs(rv - qv_)
        report(f"closed form R_{nm_} vs quadrature (200 dps)", e < mp.mpf('1e-180'), f"|R-{nm_}|={mp.nstr(e,3)}")

# ---------------- [5] exact certificates ----------------
for tgt in ("D2", "D1"):
    r = subprocess.run([PY, os.path.join(HERE, "verify_wt4_cert.py"), tgt],
                       capture_output=True, text=True, cwd=HERE, timeout=1200)
    ok = (r.returncode == 0) and (f"VERDICT: {tgt} = R PROVED" in r.stdout)
    tail = [l for l in r.stdout.splitlines() if "VERDICT" in l or "residual_nonzero" in l]
    report(f"exact certificate {tgt}", ok, "; ".join(tail[-3:]))

# ---------------- [6] k=3 offset ----------------
with mp.workdps(120):
    rr2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+rr2)
    qD1 = quad_value({k: deltas[k]["D1"] for k in ["2,1","3,0"]})
    qD2 = quad_value({k: deltas[k]["D2"] for k in ["2,1","3,0"]})
    qW4pp = quad_value(kernels["4"])
    qW43 = quad_value(kernels["3"])
    import w4_symbol as ws
    def atom_eval(a):
        if a[0] == 'li':
            z = a[2]; return mp.polylog(a[1], mp.mpf(str(z.p)) + rr2*mp.mpf(str(z.q)))
        if a[0] == 'log':
            j2, m = ws.logvec(a[1]); return mp.mpf(float(j2))*lam + mp.mpf(float(m))*L
        if a == ('pi2',): return mp.pi**2
        if a == ('z3',): return mp.zeta(3)
    thm = mp.mpf(0)
    for c, atoms in ws.conjecture_terms():
        cm = mp.mpf(str(c.p)) + rr2*mp.mpf(str(c.q))
        v = cm
        for a in atoms: v *= atom_eval(a)
        thm += v
    e3 = abs(qW43 - thm)
    delta = (qW4pp - qD1 + qD2) - thm
    report("W_4(3) quadrature == thm:W4 closed form (120 dps)", e3 < mp.mpf('1e-100'), f"|diff|={mp.nstr(e3,3)}")
    report("k=3 offset delta = P(3)-W_4(3) is nonzero", abs(delta) > mp.mpf('0.1'), f"delta={mp.nstr(delta,25)}")

print(f"\nverify_wt4.py: {sum(1 for _, ok, _ in RESULTS if ok)}/{len(RESULTS)} checks PASS  ({time.time()-t0:.0f}s)")
sys.exit(0 if all(ok for _, ok, _ in RESULTS) else 1)
'''
p=os.path.join(work,'verify_wt4.py')
open(p,'w').write(script)
print('verify_wt4.py written', len(script), 'chars')
