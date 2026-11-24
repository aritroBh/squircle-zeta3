"""WT5 step E: exact (3,2) Lie pairings of E1, E2, E3 (w5pp_symbol machinery).
Controls: W5 -> (84, -38 sqrt2); W5'' -> (15, -13 sqrt2). Saves pairings.json."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w4_symbol import *
import w5_symbol as W5

ZEROQ, ONEQ = Q(0), Q(1)
POLE = {"0": ZEROQ, "a": A, "-a": -A, "b": B, "-b": -B}
KEYS = ["4,0", "3,1", "2,2"]

def g_expression_from_terms(terms, i, j):
    factors = (LT,) * i + (MT,) * j
    out = {}
    for name, (p, q) in terms:
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            key = ((z,) + w, k)
            out[key] = out.get(key, ZEROQ) + coefficient * d
    return clean(out)

def expression_from_omegarec(rec):
    out = {}
    for key in KEYS:
        i, j = map(int, key.split(","))
        for k2, c in g_expression_from_terms(rec[key], i, j).items():
            out[k2] = out.get(k2, ZEROQ) + c
    return clean(out)

def w5pp_expression():
    om = json.loads((HERE/'omega5.json').read_text())
    return expression_from_omegarec({k: om[k] for k in KEYS})

def e_expression(which):
    dd = json.loads((HERE/'deltas_deg4.json').read_text())
    return expression_from_omegarec({k: dd[k][which] for k in KEYS})

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
        v = A_lie.get(w[:la], ZEROQ) * B_lie.get(w[la:], ZEROQ) \
            - B_lie.get(w[:lb], ZEROQ) * A_lie.get(w[lb:], ZEROQ)
        if v:
            tot += c * v
    return tot

def pairings(expr):
    S = symbol(expr)
    return [pair_ab(S, W5.L3[0], W5.L2), pair_ab(S, W5.L3[1], W5.L2)]

if __name__ == '__main__':
    ban_float(); self_check()
    # CONTROL 1: W5 pairings = (84, -38 sqrt2)
    s5 = W5.symbol()
    o5 = [pair_ab(s5, W5.L3[0], W5.L2), pair_ab(s5, W5.L3[1], W5.L2)]
    print('CONTROL W5 (3,2) pairings:', [str(v) for v in o5], '(expect 84 and -38*sqrt2)', flush=True)
    assert o5[0] == Q(84) and o5[1] == Q(0, -38), o5
    # CONTROL 2: W5'' pairings = (15, -13 sqrt2)
    p5 = pairings(w5pp_expression())
    print("CONTROL W5'' (3,2) pairings:", [str(v) for v in p5], '(expect 15 and -13*sqrt2)', flush=True)
    assert p5[0] == Q(15) and p5[1] == Q(0, -13), p5
    data = dict(control_W5=[str(v) for v in o5], control_W5pp=[str(v) for v in p5])
    for name, which in [("E1", "D1"), ("E2", "D2"), ("E3", "D3")]:
        expr = e_expression(which)
        p = pairings(expr)
        data[name] = dict(gpl_terms=len(expr), pairings=[str(v) for v in p],
                          pairings_pq=[[str(v.p), str(v.q)] for v in p])
        print(f'{name} (3,2) pairings ([X,P]xP, [Y,P]xP):', [str(v) for v in p],
              '  classical-at-symbol-level:', p[0] == ZEROQ and p[1] == ZEROQ, flush=True)
    (HERE/'pairings.json').write_text(json.dumps(data, indent=2)+'\n')
    with open(HERE/'progress.log', 'a') as f:
        f.write('step_E pairings: controls (84,-38s2),(15,-13s2) PASS; ' +
                '; '.join(f"{nm}={data[nm]['pairings']}" for nm in ['E1','E2','E3']) + '\n')
    print('saved pairings.json', flush=True)
