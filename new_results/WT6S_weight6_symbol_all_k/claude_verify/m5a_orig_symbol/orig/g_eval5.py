"""Numeric evaluation of the W6 G-expression gpl_w6.tsv vs quadrature (adapted from M4a g_eval.py)."""
import sys, json, csv, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from functools import cache
from fractions import Fraction as F
from flint import acb, arb, acb_series, acb_poly, ctx
import mpmath as mp
from qsqrt2 import QS2 as Q, A, B

ctx.prec = 800
ctx.cap = 500
N = ctx.cap
mp.mp.dps = 120

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
    n_reg = 0
    for w, k, c in load_tsv(path):
        if w and (w[0] == A or w[0] == B): n_reg += 1
        tot += qa(c)*(log2**k)*(real_g(w, ONE, A) - real_g(w, ONE, B))
    return tot, n_reg

t0 = time.time()
val, n_reg = evaluate("gpl_w6.tsv")
q = json.load(open("quad5_100.json"))
expected = mp.mpf(q["W6"])
got = mp.mpf(val.real.mid().str(110, radius=False, more=True))
err = abs(got - expected)
print(f"gpl_w6.tsv: ({time.time()-t0:.0f}s, {n_reg} endpoint-regularized words)")
print(f"   G-value = {mp.nstr(got, 60)}")
print(f"   quad    = {mp.nstr(expected, 60)}")
print(f"   |G - quad| = {mp.nstr(err, 4)}   flint rad = {val.real.rad().str(5)}", flush=True)
json.dump({"G_value": mp.nstr(got, 110), "err_vs_quad": mp.nstr(err, 6),
           "flint_rad": val.real.rad().str(8)},
          open("gpl_w6_check.json", "w"))
