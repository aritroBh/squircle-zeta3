import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""W5'' (k=5, i+j=4): G-expression, full point symbol, and the (3,2) Lie pairings of M4c.
Adapted from M5a w6_symbol.py + M4c w5_symbol.py. Control: W5 pairings must be (84, -38 sqrt2).
"""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w4_symbol import *
import w5_symbol as W5   # M4c copy: control + Lie words L2, L3

ZEROQ, ONEQ = Q(0), Q(1)
POLE = {"0": ZEROQ, "a": A, "-a": -A, "b": B, "-b": -B}

def g_expression_from_omega(omega, i, j):
    factors = (LT,) * i + (MT,) * j
    out = {}
    for name, (p, q) in omega[f"{i},{j}"]:
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            key = ((z,) + w, k)
            out[key] = out.get(key, ZEROQ) + coefficient * d
    return clean(out)

def expression_from(ijs):
    om = json.loads((HERE/'omega5.json').read_text())
    out = {}
    for ij in ijs:
        for key, c in g_expression_from_omega(om, *ij).items():
            out[key] = out.get(key, ZEROQ) + c
    return clean(out)

def w5pp_expression():
    return expression_from([(2, 2), (3, 1), (4, 0)])

def symbol(expr):
    out = {}
    for (w, k), c in expr.items():
        s = add(g_symbol(ONEQ, w, A), scale(g_symbol(ONEQ, w, B), -1))
        for _ in range(k):
            s = shuffle(s, log_symbol(Q(2)))
        out = add(out, scale(s, c))
    return out

def pair_ab(S, A_lie, B_lie):
    la = len(next(iter(A_lie))); lb = len(next(iter(B_lie)))
    tot = ZEROQ
    for w, c in S.items():
        v = A_lie.get(w[:la], ZEROQ) * B_lie.get(w[la:], ZEROQ) \\
            - B_lie.get(w[:lb], ZEROQ) * A_lie.get(w[lb:], ZEROQ)
        if v:
            tot += c * v
    return tot

if __name__ == '__main__':
    ban_float(); self_check()
    # CONTROL: W5 pairings = (84, -38 sqrt2)
    s5 = W5.symbol()
    o5 = [pair_ab(s5, W5.L3[0], W5.L2), pair_ab(s5, W5.L3[1], W5.L2)]
    print('CONTROL W5 (3,2) pairings:', [str(v) for v in o5], '(expect 84 and -38*sqrt2)', flush=True)
    assert o5[0] == Q(84) and o5[1] == Q(0, -38), o5
    expr = w5pp_expression()
    print("W5'' gpl terms:", len(expr), flush=True)
    assert all(len(w) + k == 5 for (w, k) in expr)
    S = symbol(expr)
    p = [pair_ab(S, W5.L3[0], W5.L2), pair_ab(S, W5.L3[1], W5.L2)]
    print("W5'' (3,2) pairings ([X,P]xP, [Y,P]xP):", [str(v) for v in p], flush=True)
    # save tsv
    lines = ["# Sum c*log(2)^k*(G_1(word;a)-G_1(word;b)); first letter outermost. W5'' (k=5 weight-5 part).",
             'coefficient_p\\tcoefficient_q\\tk\\tword']
    for (w, k), c in expr.items():
        lines.append(f'{c.p}\\t{c.q}\\t{k}\\t{w}')
    (HERE/'gpl_w5pp.tsv').write_text('\\n'.join(lines)+'\\n')
    data = dict(gpl_terms=len(expr), control_W5=[str(v) for v in o5],
                pairings=[str(v) for v in p],
                symbol={''.join(map(str, w)): [str(c.p), str(c.q)] for w, c in sorted(S.items())})
    (HERE/'symbol5pp.json').write_text(json.dumps(data, indent=2)+'\\n')
    print('saved gpl_w5pp.tsv and symbol5pp.json', flush=True)
'''
open(os.path.join(m5b,'w5pp_symbol.py'),'w').write(script)
# also the evaluator: copy g_eval5.py pattern for the new tsv
ev = '''"""Numeric evaluation of gpl_w5pp.tsv vs quadrature (adapted from M5a g_eval5.py)."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M5a_next_coefficient")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from g_eval5 import evaluate, mp
t0 = time.time()
val, n_reg = evaluate("gpl_w5pp.tsv")
q = json.load(open("quad5_100.json"))
expected = mp.mpf(q["W5pp"])
got = mp.mpf(val.real.mid().str(110, radius=False, more=True))
err = abs(got - expected)
print(f"gpl_w5pp.tsv: ({time.time()-t0:.0f}s, {n_reg} endpoint-regularized words)")
print(f"   G-value = {mp.nstr(got, 60)}")
print(f"   quad    = {mp.nstr(expected, 60)}")
print(f"   |G - quad| = {mp.nstr(err, 4)}   flint rad = {val.real.rad().str(5)}", flush=True)
json.dump({"G_value": mp.nstr(got, 110), "err_vs_quad": mp.nstr(err, 6),
           "flint_rad": val.real.rad().str(8)}, open("gpl_w5pp_check.json", "w"))
'''
open(os.path.join(m5b,'g_eval5pp.py'),'w').write(ev)
print('written w5pp_symbol.py and g_eval5pp.py')
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','w5pp_symbol.py'],capture_output=True,text=True,cwd=m5b,timeout=1200)
print(r.stdout[-2500:]); print(r.stderr[-800:])
