base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step B: exact G-expressions (gpl TSVs) for E1, E2, E3 from deltas_deg4.json;
control: rebuild W5pp G-expression from kernels_deg4.json k=5 and compare to M5b gpl_w5pp.tsv EXACTLY."""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from collections import defaultdict
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
from w4_symbol import log_product, LT, MT, clean

t0 = time.time()
ZERO = Q(0)
POLE = {"0": ZERO, "a": A, "-a": -A, "b": B, "-b": -B}
KEYS = ["4,0", "3,1", "2,2"]

def g_expression(i, j, terms):
    out = defaultdict(lambda: Q(0))
    factors = (LT,) * i + (MT,) * j
    for name, (p, q) in terms:
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            out[((z,) + w, k)] += coefficient * d
    return clean(out)

def combine(*dicts):
    out = defaultdict(lambda: Q(0))
    for d in dicts:
        for k, v in d.items():
            out[k] += v
    return clean(out)

def build(kerns):
    return combine(g_expression(4, 0, kerns["4,0"]), g_expression(3, 1, kerns["3,1"]), g_expression(2, 2, kerns["2,2"]))

def save_tsv(gex, path, header):
    lines = [header, "coefficient_p\\tcoefficient_q\\tk\\tword"]
    for (w, k), c in gex.items():
        lines.append(f"{c.p}\\t{c.q}\\t{k}\\t{w}")
    open(path, "w").write("\\n".join(lines) + "\\n")

kernels = json.load(open(os.path.join(HERE, "kernels_deg4.json")))
deltas = json.load(open(os.path.join(HERE, "deltas_deg4.json")))

gE1 = build({k: deltas[k]["D1"] for k in KEYS})
gE2 = build({k: deltas[k]["D2"] for k in KEYS})
gE3 = build({k: deltas[k]["D3"] for k in KEYS})
gW5 = build(kernels["5"])

assert all(len(w) + k == 5 for (w, k) in gE1)
assert all(len(w) + k == 5 for (w, k) in gE2)
assert all(len(w) + k == 5 for (w, k) in gE3)
assert all(len(w) + k == 5 for (w, k) in gW5)
print(f"E1 terms={len(gE1)}  E2 terms={len(gE2)}  E3 terms={len(gE3)}  W5pp terms={len(gW5)}  t={time.time()-t0:.1f}s")

save_tsv(gE1, os.path.join(HERE, "e1.tsv"),
         "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. E1 = W5 slope (K^(6)-K^(5)).")
save_tsv(gE2, os.path.join(HERE, "e2.tsv"),
         "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. E2 = W5 second difference (K^(7)-2K^(6)+K^(5)).")
save_tsv(gE3, os.path.join(HERE, "e3.tsv"),
         "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. E3 = W5 third difference (K^(8)-3K^(7)+3K^(6)-K^(5)).")

# CONTROL: compare rebuilt W5pp to M5b's gpl_w5pp.tsv as exact dicts
import csv
def load_tsv_dict(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\\t"))
    out = defaultdict(lambda: Q(0))
    for r in rows:
        out[(r["word"], int(r["k"]))] += Q(F(r["coefficient_p"]), F(r["coefficient_q"]))
    return dict(out)
ref = load_tsv_dict(os.path.join(HERE, "gpl_w5pp.tsv"))
mine = {(str(w), k): c for (w, k), c in gW5.items()}
assert mine == ref, (len(mine), len(ref), list(set(mine.items()) ^ set(ref.items()))[:5])
print(f"CONTROL: rebuilt W5pp G-expression == M5b gpl_w5pp.tsv EXACTLY ({len(ref)} terms)")

for nm, g in [("E1", gE1), ("E2", gE2), ("E3", gE3)]:
    al = sorted({x for (w, k) in g for x in w}, key=float)
    print(f"{nm} alphabet: {[str(x) for x in al]}")
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_B gexpr: E1 {len(gE1)} terms, E2 {len(gE2)} terms, E3 {len(gE3)} terms; control W5pp==gpl_w5pp.tsv EXACT; t={time.time()-t0:.1f}s\\n")
print("DONE", time.time()-t0, "s")
'''
open(W+'wt5_B_gexpr.py','w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_B_gexpr.py'], cwd=W, capture_output=True, text=True, timeout=900)
print(r.stdout); print('STDERR:', r.stderr[-1500:])
