import os
work='<project>/research/08_theory/WT4_weight4_all_k'
script = r'''"""WT4 step C: d1.tsv / d2.tsv G-values (flint Taylor evaluator with endpoint regularization,
ctx.prec=700 bits) vs direct mpmath quadrature of the kernel integrals at dps 130.
CONTROL: same code on kernels k=4 reproduces gpl_w4pp.tsv / quad W4'' (M4a quad_90.json)."""
import sys, os, json, csv, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from functools import cache
from fractions import Fraction as F
from flint import acb, arb, acb_series, acb_poly, ctx
import mpmath as mp
from qsqrt2 import QS2 as Q, A, B

ctx.prec = 700
ctx.cap = 420
N = ctx.cap
mp.mp.dps = 130
print("mp.mp.dps =", mp.mp.dps, " ctx.prec(bits) =", ctx.prec, flush=True)

def qa(c):
    return acb(str(c.p)) + acb(2).sqrt()*acb(str(c.q))

@cache
def shd(u, v):
    if not u: return ((v, 1),)
    if not v: return ((u, 1),)
    d = {}
    for w, m in shd(u[1:], v):
        d[(u[0],)+w] = d.get((u[0],)+w, 0) + m
    for w, m in shd(u, v[1:]):
        d[(v[0],)+w] = d.get((v[0],)+w, 0) + m
    return tuple(d.items())

@cache
def reg_words(w, e):
    if not w or w[0] != e:
        return ((w, F(1)),)
    terms = dict(shd((e,), w[1:]))
    cw = terms.pop(w)
    d = {}
    for u, m in terms.items():
        for t, cf in reg_words(u, e):
            d[t] = d.get(t, F(0)) + F(-m, cw)*cf
    return tuple((t, cf) for t, cf in d.items() if cf)

@cache
def real_series(poles):
    if not poles: return acb_series([1])
    assert poles[-1] != Q(0)
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
    value = taylor(w, c, e)
    if value is not None: return value
    value = taylor(tuple(reversed(w)), e, c)
    if value is not None: return (-1)**len(w)*value
    mid = (c+e)/2
    return sum((real_g(w[:k], mid, e)*real_g(w[k:], c, mid) for k in range(len(w)+1)), acb(0))

def real_g(w, c, e):
    tot = acb(0)
    for u, cf in reg_words(w, e):
        tot += acb(cf.numerator)/cf.denominator*real_g_raw(u, c, e)
    return tot

def parse_word(s):
    vals = {"0": Q(0), "(1+1*sqrt2)": A, "(-1-1*sqrt2)": -A,
            "(-1+1*sqrt2)": B, "(1-1*sqrt2)": -B}
    s = s.strip(); assert s.startswith("(") and s.endswith(")")
    body = s[1:-1].strip()
    if not body: return ()
    parts, depth, cur = [], 0, ""
    for ch in body:
        if ch == "(": depth += 1
        if ch == ")": depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur.strip()); cur = ""
        else:
            cur += ch
    if cur.strip(): parts.append(cur.strip())
    return tuple(vals[p] for p in parts)

def load_tsv(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    return [(parse_word(r["word"]), int(r["k"]), Q(F(r["coefficient_p"]), F(r["coefficient_q"]))) for r in rows]

def evaluate(path):
    ONE = Q(1); log2 = acb(2).log(); tot = acb(0)
    for w, k, c in load_tsv(path):
        tot += qa(c)*(log2**k)*(real_g(w, ONE, A) - real_g(w, ONE, B))
    return tot

# ------------- mpmath quadrature of kernel integrals -------------
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
def quad_value(kerns):
    tot = mp.mpf(0); err = mp.mpf(0)
    for key, i, j in [("2,1", 2, 1), ("3,0", 3, 0)]:
        ts = [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in kerns[key]]
        def f(v, ts=ts, i=i, j=j):
            l1 = mp.log((v - bN)*(v + aN)/(4*v))
            l2 = mp.log((aN - v)*(v + bN)/(4*v))
            return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
        val, e = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
        tot += val; err = max(err, e)
    return tot, err

deltas = json.load(open(os.path.join(HERE, "deltas_deg3.json")))
kernels = json.load(open(os.path.join(HERE, "kernels_deg3.json")))
targets = [
    ("D1", "d1.tsv", {k: deltas[k]["D1"] for k in ["2,1", "3,0"]}),
    ("D2", "d2.tsv", {k: deltas[k]["D2"] for k in ["2,1", "3,0"]}),
    ("W4pp(control)", os.path.join("<project>/research/08_theory/M4a_decomposition", "gpl_w4pp.tsv"), kernels["4"]),
]
res = {}
for name, tsv, kerns in targets:
    t1 = time.time()
    gv = evaluate(tsv)
    qv, qe = quad_value(kerns)
    got = mp.mpf(gv.real.mid().str(150, radius=False, more=True))
    diff = abs(got - qv)
    rad = gv.real.rad()
    print(f"{name}: ({time.time()-t1:.0f}s)")
    print(f"   G-value = {mp.nstr(got, 70)}")
    print(f"   quad    = {mp.nstr(qv, 70)}  (quad err_est {mp.nstr(qe,3)})")
    print(f"   |G - quad| = {mp.nstr(diff, 4)}   flint rad = {rad.str(5)}   at mp.mp.dps={mp.mp.dps}", flush=True)
    res[name] = {"G_value": mp.nstr(got, 130), "quad": mp.nstr(qv, 130), "diff": mp.nstr(diff, 6),
                 "flint_rad": rad.str(8), "dps": mp.mp.dps}
# control vs M4a quad_90
q90 = json.load(open("<project>/research/08_theory/M4a_decomposition/quad_90.json"))
w4pp_ref = mp.mpf(q90["W4pp"])
ctl = abs(mp.mpf(res["W4pp(control)"]["quad"]) - w4pp_ref)
print(f"CONTROL quad(W4pp) - M4a quad_90 W4pp = {mp.nstr(ctl, 4)} at mp.mp.dps={mp.mp.dps}")
res["control_quad_vs_m4a"] = mp.nstr(ctl, 6)
json.dump(res, open(os.path.join(HERE, "check_tsv_vs_quad.json"), "w"), indent=1)
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_C check: D1 |G-quad|={res['D1']['diff']}, D2 |G-quad|={res['D2']['diff']}, W4pp ctl |G-quad|={res['W4pp(control)']['diff']}, ctl-vs-m4a={res['control_quad_vs_m4a']} at dps {mp.mp.dps}\n")
print("DONE")
'''
p=os.path.join(work,'wt4_C_check.py')
open(p,'w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','19',PY,'wt4_C_check.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR:', r.stderr[-1500:])
