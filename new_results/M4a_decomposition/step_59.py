code = '''"""G-expressions for W5 and W4'' (k=4), generalising w4_symbol.w4_g_expression.
Input: omega4.json (exact pole/residue parse of ibp_reduce k=4 output).
Output: {(z,)+w, k): c} dicts + TSV files in gpl_expression.tsv format."""
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
    """terms: list of (pole_name, (p_str, q_str)); piece int l1^i l2^j * sum c/(v-z)."""
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
    om = json.load(open("omega4.json"))
    return {tuple(map(int, k.split(","))): v for k, v in om.items()}

def w5_g_expression():
    om = load()
    return combine(g_expression(2, 2, om[(2, 2)]),
                   g_expression(3, 1, om[(3, 1)]),
                   g_expression(4, 0, om[(4, 0)]))

def w4pp_g_expression():
    om = load()
    return combine(g_expression(2, 1, om[(2, 1)]),
                   g_expression(3, 0, om[(3, 0)]))

def save_tsv(gex, path, header):
    lines = [header, "coefficient_p\tcoefficient_q\tk\tword"]
    for (w, k), c in gex.items():
        lines.append(f"{c.p}\t{c.q}\t{k}\t{w}")
    open(path, "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    w5 = w5_g_expression()
    w4 = w4pp_g_expression()
    save_tsv(w5, "gpl_w5.tsv", "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. W5 (k=4 weight-5 part).")
    save_tsv(w4, "gpl_w4pp.tsv", "# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. W4'' (k=4 weight-4 part).")
    print("W5 terms:", len(w5), " max word length:", max(len(w) for (w, k) in w5))
    print("W4'' terms:", len(w4), " max word length:", max(len(w) for (w, k) in w4))
    # sanity: every entry has len(w)+k == weight
    assert all(len(w) + k == 5 for (w, k) in w5)
    assert all(len(w) + k == 4 for (w, k) in w4)
    print("weight check: all entries satisfy len(word)+k == weight: OK")
'''
open("gexpr.py", "w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable, "gexpr.py"], capture_output=True, text=True, timeout=1200)
print(r.stdout)
print("ERR", r.stderr[-800:])
