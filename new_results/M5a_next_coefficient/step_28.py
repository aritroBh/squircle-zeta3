work = '<project>/research/08_theory/M5a_next_coefficient'
open(f'{work}/gexpr5.py', 'w').write('''"""G-expression for W6 (k=5, weight-6 part), generalising M4a gexpr.py."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from collections import defaultdict
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
from w4_symbol import log_product, LT, MT, clean

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

def load():
    om = json.load(open("omega5.json"))
    return {tuple(map(int, k.split(","))): v for k, v in om.items()}

def w6_g_expression():
    om = load()
    return combine(g_expression(3, 2, om[(3, 2)]),
                   g_expression(4, 1, om[(4, 1)]),
                   g_expression(5, 0, om[(5, 0)]))

def save_tsv(gex, path, header):
    lines = [header, "coefficient_p\\tcoefficient_q\\tk\\tword"]
    for (w, k), c in gex.items():
        lines.append(f"{c.p}\\t{c.q}\\t{k}\\t{w}")
    open(path, "w").write("\\n".join(lines) + "\\n")

if __name__ == "__main__":
    w6 = w6_g_expression()
    save_tsv(w6, "gpl_w6.tsv", "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. W6 (k=5 weight-6 part).")
    print("W6 terms:", len(w6), " max word length:", max(len(w) for (w, k) in w6))
    assert all(len(w) + k == 6 for (w, k) in w6)
    print("weight check len(word)+k == 6: OK")
    al6 = sorted({x for (w, k) in w6 for x in w}, key=float)
    print("W6 alphabet:", [str(x) for x in al6])
    firstletters = sorted({w[0] for (w, k) in w6 if w}, key=float)
    print("first letters:", [str(x) for x in firstletters])
''')
import subprocess, time
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'gexpr5.py'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
if 'alphabet' in r.stdout:
    with open(f'{work}/progress.log', 'a') as f:
        f.write("step: gexpr5.py built gpl_w6.tsv (log_product G-expression of W6)\n")