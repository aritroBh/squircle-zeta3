base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step L: exact cubic law at the G-expression level: for k=6..14, the W5(k) G-expression
rebuilt from kernels_deg4.json[k] equals W5pp + (k-5)E1 + C(k-5,2)E2 + C(k-5,3)E3 as exact TSV dicts."""
import sys, os, json, csv, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from collections import defaultdict
from fractions import Fraction as F
from math import comb
from qsqrt2 import QS2 as Q
from w4_symbol import log_product, LT, MT, clean

t0 = time.time()
ZERO = Q(0)
POLE = {"0": ZERO, "a": A, "-a": -A, "b": B, "-b": -B}
KEYS = ["4,0", "3,1", "2,2"]
from qsqrt2 import A, B

def g_expression(i, j, terms):
    out = defaultdict(lambda: Q(0))
    factors = (LT,) * i + (MT,) * j
    for name, (p, q) in terms:
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            out[((z,) + w, k)] += coefficient * d
    return clean(out)

def build(kerns):
    out = defaultdict(lambda: Q(0))
    for key in KEYS:
        i, j = map(int, key.split(","))
        for kk, c in g_expression(i, j, kerns[key]).items():
            out[kk] += c
    return clean(out)

def load_tsv_dict(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\\t"))
    out = defaultdict(lambda: Q(0))
    for r in rows:
        out[(r["word"], int(r["k"]))] += Q(F(r["coefficient_p"]), F(r["coefficient_q"]))
    return dict(out)

kernels = json.load(open(os.path.join(HERE, "kernels_deg4.json")))
g0 = {(str(w), k): c for (w, k), c in build(kernels["5"]).items()}
gE = {1: load_tsv_dict(os.path.join(HERE, "e1.tsv")),
      2: load_tsv_dict(os.path.join(HERE, "e2.tsv")),
      3: load_tsv_dict(os.path.join(HERE, "e3.tsv"))}
ok = True
for k in range(6, 15):
    gk = {(str(w), k2): c for (w, k2), c in build(kernels[str(k)]).items()}
    pred = defaultdict(lambda: Q(0))
    for key in set(g0) | set(gE[1]) | set(gE[2]) | set(gE[3]):
        pred[key] = (g0.get(key, ZERO) + (k-5)*gE[1].get(key, ZERO)
                     + comb(k-5, 2)*gE[2].get(key, ZERO) + comb(k-5, 3)*gE[3].get(key, ZERO))
    pred = {kk: c for kk, c in pred.items() if c != ZERO}
    same = pred == gk
    ok = ok and same
    print(f"k={k}: G-expr terms={len(gk)} cubic-law-exact={same}  t={time.time()-t0:.1f}s", flush=True)
assert ok
print("CUBIC G-EXPRESSION LAW exact for k=6..14: PASS")
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_L: cubic G-expression law exact k=6..14 (TSV-level); t={time.time()-t0:.1f}s\\n")
'''
open(W+'wt5_L_cubicg.py','w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_L_cubicg.py'], cwd=W, capture_output=True, text=True, timeout=1200)
print(r.stdout); print('STDERR:', r.stderr[-500:])
