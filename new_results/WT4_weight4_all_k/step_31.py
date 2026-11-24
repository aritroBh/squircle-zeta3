import os
work='<project>/research/08_theory/WT4_weight4_all_k'
script = r'''"""WT4 step B: exact G-expressions (gpl TSVs) for D1, D2 from deltas_deg3.json;
control: rebuild W4pp G-expression from kernels_deg3.json k=4 and compare to M4a gpl_w4pp.tsv EXACTLY."""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from collections import defaultdict
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
from w4_symbol import log_product, LT, MT, clean

t0 = time.time()
ZERO = Q(0)
POLE = {"0": ZERO, "a": A, "-a": -A, "b": B, "-b": -B}

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
    return combine(g_expression(2, 1, kerns["2,1"]), g_expression(3, 0, kerns["3,0"]))

def save_tsv(gex, path, header):
    lines = [header, "coefficient_p\tcoefficient_q\tk\tword"]
    for (w, k), c in gex.items():
        lines.append(f"{c.p}\t{c.q}\t{k}\t{w}")
    open(path, "w").write("\n".join(lines) + "\n")

kernels = json.load(open(os.path.join(HERE, "kernels_deg3.json")))
deltas = json.load(open(os.path.join(HERE, "deltas_deg3.json")))

gD1 = build(deltas and {k: deltas[k]["D1"] for k in ["2,1", "3,0"]})
gD2 = build({k: deltas[k]["D2"] for k in ["2,1", "3,0"]})
gW4 = build(kernels["4"])

assert all(len(w) + k == 4 for (w, k) in gD1)
assert all(len(w) + k == 4 for (w, k) in gD2)
assert all(len(w) + k == 4 for (w, k) in gW4)
print(f"D1 terms={len(gD1)}  D2 terms={len(gD2)}  W4pp terms={len(gW4)}  t={time.time()-t0:.1f}s")

save_tsv(gD1, os.path.join(HERE, "d1.tsv"),
         "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. D1 = W4 slope (K^(5)-K^(4)).")
save_tsv(gD2, os.path.join(HERE, "d2.tsv"),
         "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. D2 = W4 second difference (K^(6)-2K^(5)+K^(4)).")

# CONTROL: compare rebuilt W4pp to M4a's gpl_w4pp.tsv as exact dicts
import csv
def load_tsv_dict(path):
    rows = list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))
    out = defaultdict(lambda: Q(0))
    for r in rows:
        out[(r["word"], int(r["k"]))] += Q(F(r["coefficient_p"]), F(r["coefficient_q"]))
    return dict(out)
ref = load_tsv_dict("<project>/research/08_theory/M4a_decomposition/gpl_w4pp.tsv")
mine = {(str(w), k): c for (w, k), c in gW4.items()}
assert mine == ref, (set(mine.items()) ^ set(ref.items()))
print(f"CONTROL: rebuilt W4pp G-expression == M4a gpl_w4pp.tsv EXACTLY ({len(ref)} terms)")

# alphabets
for nm, g in [("D1", gD1), ("D2", gD2)]:
    al = sorted({x for (w, k) in g for x in w}, key=float)
    print(f"{nm} alphabet: {[str(x) for x in al]}")
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write(f"step_B gexpr: D1 {len(gD1)} terms, D2 {len(gD2)} terms; control W4pp==gpl_w4pp.tsv EXACT; t={time.time()-t0:.1f}s\n")
print("DONE", time.time()-t0, "s")
'''
p=os.path.join(work,'wt4_B_gexpr.py')
open(p,'w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','19',PY,'wt4_B_gexpr.py'], cwd=work, capture_output=True, text=True, timeout=900)
print(r.stdout)
print('STDERR:', r.stderr[-2000:])
