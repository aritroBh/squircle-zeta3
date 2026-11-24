verify='''"""M5c verify W6 = R6 + N6 (tag='' or w6x..):
1. EXACT symbol of R6 equals symbol(W6) - symbol(N6).
2. EXACT (3,3) and (4,2) pairings of R6 + N6 - W6 are 0.
3. Value of R6 + N6 - W6 at >=3000 digits vs w6_11800 (independent quadrature).
Usage: -B w6c_verify.py [tag]
"""
import json, math, sys, pickle
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp
from w4_symbol import Q, R, A, B, ban_float, g_symbol, add, scale, concat
import w5pp_symbol as WS
from m5b_common import ZERO, ONE
from w6_solve import atoms_w6, n6_terms, mpl_word, target_symbol6, parse_qs2, BASE_LI6, pool_zs
from w6_values import leaf_recipes

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

def li6_args_for(tag):
    if not tag:
        return list(BASE_LI6)
    idx = [int(i) for i in tag[3:].split('_')]
    pool = pool_zs()
    li6 = list(BASE_LI6)
    for i in idx:
        if pool[i] not in li6:
            li6.append(pool[i])
    return li6

def main():
    ban_float()
    tag = sys.argv[1] if len(sys.argv) > 1 else ''
    d = json.loads((HERE/'w6_closed_form.json').read_text())
    names = d['names']; ca = [F(c) for c in d['ca']]; cb = [F(c) for c in d['cb']]
    basis = atoms_w6(li6_args_for(tag))
    symmap = {nm: s for nm, s, _ in basis}
    t = {}
    for nm, a, b in zip(names, ca, cb):
        t = add(t, scale(symmap[nm], Q(a, b)))
    ts = target_symbol6()
    assert t == ts, 'symbol mismatch'
    print('1. EXACT symbol of R6 equals symbol(W6) - symbol(N6): PASS', flush=True)
    s_w6 = WS.symbol(WS.expression_from([(3, 2), (4, 1), (5, 0)]))
    s_n6 = {}
    for m, n, x, y, c in n6_terms():
        s_n6 = add(s_n6, scale(g_symbol(ZERO, mpl_word(m, n, x, y), ONE), c))
    diff = add(t, s_n6, scale(s_w6, -1))
    assert not diff, 'R6+N6-W6 symbol nonzero'
    p33 = pair33(diff); p42 = pair42(diff)
    assert p33 == ZERO and all(v == ZERO for v in p42)
    print('2. EXACT (3,3) and (4,2) pairings of R6+N6-W6 = 0: PASS', flush=True)
    # 3. value at >=3000 digits (3500-digit caches)
    mp.mp.dps = 3500 + 120
    suf = '' if not tag else '_' + tag.replace('w6x','x')
    recs = leaf_recipes()
    saved = pickle.loads((HERE/f'w6_leaves_3500{ suf if False else ""}.pkl').read_bytes())
    if isinstance(saved, dict):
        saved = list(saved.values())
    leaves = {r: mp.mpf(v) for r, v in zip(recs, saved)}
    if tag:
        extra = pickle.loads((HERE/f'w6_leaves_extra_{tag}.pkl').read_bytes())
        for rec, v in extra.items():
            leaves[rec] = mp.mpf(v)
    vals = [mp.mpf(v) for v in pickle.loads((HERE/f'w6_atomvals_3500{suf}.pkl').read_bytes())]
    name2val = {b[0]: v for b, v in zip(basis, vals)}
    r2 = mp.sqrt(2)
    total = mp.mpf(0)
    for nm, a, b in zip(names, ca, cb):
        total += (mp.mpf(a.numerator)/a.denominator + r2*mp.mpf(b.numerator)/b.denominator)*name2val[nm]
    n6v = mp.mpf(0)
    for m, n, x, y, c in n6_terms():
        n6v += (mp.mpf(c.p.numerator)/c.p.denominator + r2*mp.mpf(c.q.numerator)/c.q.denominator)*leaves[('mpl', m, n, x, y)]
    total += n6v
    w6 = mp.mpf(json.loads((HERE/'w6_11800.json').read_text())['value_full'])
    res = total - w6
    print(f'3. R6 + N6 - W6(11800 bits): {mp.nstr(res, 8)}  (expect ~1e-3400 or better)', flush=True)
    print('ALL W6 VERIFICATION CHECKS DONE')

main()
'''
open('w6c_verify.py','w').write(verify)
import py_compile; py_compile.compile('w6c_verify.py', doraise=True)
print('w6c_verify.py ok')
