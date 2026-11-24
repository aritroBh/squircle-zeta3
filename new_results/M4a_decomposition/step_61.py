import subprocess, sys, datetime
code = r'''"""Numeric evaluation of real-letter G-expressions (adapted from A3c check_conversion.py real_g),
self-contained (no writes into reference dirs). Check gpl_w5.tsv / gpl_w4pp.tsv vs quadrature."""
import sys, json, csv
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from functools import cache
from fractions import Fraction as F
from flint import acb, arb, acb_series, acb_poly, ctx
import mpmath as mp
from qsqrt2 import QS2 as Q, A, B

ctx.prec = 700
ctx.cap = 420
N = ctx.cap
mp.mp.dps = 100

def qa(c):
    return acb(str(c.p)) + acb(2).sqrt()*acb(str(c.q))

@cache
def real_series(poles):
    if not poles: return acb_series([1])
    assert poles[-1] != Q(0)
    tail = real_series(poles[1:])
    if poles[0] == Q(0): return acb_series([0]+[tail[j]/j for j in range(1, N)])
    return (tail/acb_series([-qa(poles[0]), 1])).integral()

@cache
def real_g(w, c, e):
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

def load_tsv(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    return [((eval(row["word"], {}, {"__builtins__": {}})), int(row["k"]),
             Q(F(row["coefficient_p"]), F(row["coefficient_q"]))) for row in rows]

def parse_qs2(s):
    return s

import ast
def parse_word(s):
    # word is repr of tuple of QS2; rebuild via mapping of the printed atoms
    import re
    vals = {"0": Q(0), "(1+1*sqrt2)": A, "(-1-1*sqrt2)": -A,
            "(-1+1*sqrt2)": B, "(1-1*sqrt2)": -B}
    s = s.strip()
    assert s.startswith("(") and s.endswith(")")
    body = s[1:-1].strip()
    if not body: return ()
    # split top-level commas: atoms themselves are parenthesised or '0'
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

def load_tsv2(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    out = []
    for row in rows:
        w = parse_word(row["word"])
        out.append((w, int(row["k"]), Q(F(row["coefficient_p"]), F(row["coefficient_q"]))))
    return out

def evaluate(path):
    terms = load_tsv2(path)
    ONE = Q(1)
    log2 = acb(2).log()
    tot = acb(0)
    for w, k, c in terms:
        g = real_g(w, ONE, A) - real_g(w, ONE, B)
        tot += qa(c) * (log2**k) * g
    return tot

for path, quadkey in [("gpl_w4pp.tsv", "W4pp"), ("gpl_w5.tsv", "W5")]:
    val = evaluate(path)
    q = json.load(open("quad_90.json"))
    expected = mp.mpf(q[quadkey])
    got = mp.mpf(val.real.mid().str(90, radius=False, more=True))
    err = abs(got - expected)
    rad = val.real.rad()
    print(f"{path}: G-value = {mp.nstr(got, 60)}")
    print(f"   quad    = {mp.nstr(expected, 60)}")
    print(f"   |G - quad| = {mp.nstr(err, 4)}   flint rad = {rad.str(5)}", flush=True)
'''
open("g_eval.py", "w").write(code)
r = subprocess.run([sys.executable, "g_eval.py"], capture_output=True, text=True, timeout=2700)
print(r.stdout)
print("ERR", r.stderr[-800:])
