"""Verify W6 = R6 + N6:
1. EXACT symbol of R6 equals symbol(W6) - symbol(N6).
2. EXACT (3,3) and (4,2) pairings of R6 + N6 - W6 are 0.
3. Value of R6 + N6 - W6 at >=1600 digits vs w6_5800 (independent quadrature).
"""
import json, math, sys, itertools
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from flint import ctx
import mpmath as mp
from w4_symbol import Q, R, A, B, ban_float, g_symbol, add, scale
from w4_symbol import concat
import w5_symbol as W
import w5pp_symbol as WS
from m5b_common import ZERO, ONE, li32_word
from w6_solve import atoms_w6, n6_terms, mpl_word, BASE_LI6, target_symbol6
from w6_lll import build, set_prec
from w6_values import leaf_recipes
import pickle

def bracket(x, y):
    return add(concat(x, y), scale(concat(y, x), -1))

Xl, Yl = {(0,): ONE}, {(1,): ONE}
Pl = bracket(Xl, Yl)
B3 = [bracket(Xl, Pl), bracket(Yl, Pl)]
B4 = [bracket(Xl, B3[0]), bracket(Yl, B3[0]), bracket(Yl, B3[1])]

def pair_ab(S, A_lie, B_lie):
    la = len(next(iter(A_lie))); lb = len(next(iter(B_lie)))
    tot = ZERO
    for w, c in S.items():
        v = A_lie.get(w[:la], ZERO)*B_lie.get(w[la:], ZERO) - B_lie.get(w[:lb], ZERO)*A_lie.get(w[lb:], ZERO)
        if v:
            tot += c*v
    return tot

def pair33(S):
    return pair_ab(S, B3[0], B3[1])

def pair42(S):
    return [pair_ab(S, A4, Pl) for A4 in B4]

def main():
    ban_float()
    d = json.loads((HERE/'w6_closed_form.json').read_text())
    names = d['names']; ca = [F(c) for c in d['ca']]; cb = [F(c) for c in d['cb']]
    basis = atoms_w6(BASE_LI6)
    symmap = {nm: s for nm, s, _ in basis}
    recmap = {nm: rec for nm, _, rec in basis}
    # 1. exact symbol
    t = {}
    for nm, a, b in zip(names, ca, cb):
        t = add(t, scale(symmap[nm], Q(a, b)))
    ts = target_symbol6()
    assert t == ts, 'symbol mismatch'
    print('1. EXACT symbol of R6 equals symbol(W6) - symbol(N6): PASS', flush=True)
    # 2. pairings of R6 + N6 - W6 = 0 (their symbols are equal; compute directly)
    s_w6 = WS.symbol(WS.expression_from([(3, 2), (4, 1), (5, 0)]))
    s_n6 = {}
    for m, n, x, y, c in n6_terms():
        s_n6 = add(s_n6, scale(g_symbol(ZERO, mpl_word(m, n, x, y), ONE), c))
    s_r6 = t
    diff = add(s_r6, s_n6, scale(s_w6, -1))
    assert not diff, 'R6+N6-W6 symbol nonzero'
    p33 = pair33(diff); p42 = pair42(diff)
    assert p33 == ZERO and all(v == ZERO for v in p42)
    print('2. EXACT (3,3) and (4,2) pairings of R6+N6-W6 = 0: PASS', flush=True)
    # 3. value at >=1600 digits
    set_prec(1600)
    recs = leaf_recipes()
    saved = pickle.loads((HERE/'w6_leaves_1650.pkl').read_bytes())
    if isinstance(saved, dict):
        saved = list(saved.values())
    leaves = {r: mp.mpf(v) for r, v in zip(recs, saved)}
    vals = [mp.mpf(v) for v in pickle.loads((HERE/'w6_atomvals_1650.pkl').read_bytes())]
    name2val = {b[0]: v for b, v in zip(basis, vals)}
    total = mp.mpf(0)
    r2 = mp.sqrt(2)
    for nm, a, b in zip(names, ca, cb):
        total += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*name2val[nm]
    n6v = mp.mpf(0)
    for m, n, x, y, c in n6_terms():
        n6v += (mp.mpf(c.p.numerator)/c.p.denominator + r2*mp.mpf(c.q.numerator)/c.q.denominator)*leaves[('mpl', m, n, x, y)]
    total += n6v
    w6 = mp.mpf(json.loads((HERE/'w6_5800.json').read_text())['value_full'])
    res = total - w6
    print(f'3. R6 + N6 - W6(5800 bits): {mp.nstr(res, 8)}  (expect ~1e-1600)', flush=True)
    print('ALL W6 VERIFICATION CHECKS DONE')

main()
