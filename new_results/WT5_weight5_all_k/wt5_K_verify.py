"""WT5 step K: verify E_i = R_i + N_i for e1, e2, e3:
1. EXACT full point symbol of R_i equals symbol(E_i) - symbol(N_i).
2. EXACT (3,2) cobracket: each Li5 -> 0; shuffle-product control; reduce(N_i - target(E_i)) = 0.
3. Value at ~1400 digits vs e{i}_5800 quadrature (1715 digits).
Usage: -B wt5_K_verify.py <e1|e2|e3> [closed_form_suffix]
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, acb, ctx
import mpmath as mp
from w4_symbol import Q, R, A, B, ban_float, g_symbol, add, scale
import w5_symbol as W
import wt5_E_pairings as WS
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word
from w4_right_coproduct import P, PE
from w5_search import evaluate
from m5b_common import ZERO, ONE, pool_zs, atoms_w5, li32_word, qnum

def main():
    name = sys.argv[1]
    suffix = sys.argv[2] if len(sys.argv) > 2 else ''
    which = {'e1': 'D1', 'e2': 'D2', 'e3': 'D3'}[name]
    ban_float()
    d = json.loads((HERE/f'{name}_closed_form{suffix}.json').read_text())
    names = d['names']; ca = [F(c) for c in d['ca']]; cb = [F(c) for c in d['cb']]
    aL, bL = F(d['aL']), F(d['bL'])
    assert aL == 0 and bL == 0, 'L(chi8,5) present'
    zs = pool_zs()
    basis = atoms_w5(zs)
    symmap = {nm: s for nm, s, _ in basis}
    recmap = {nm: rec for nm, _, rec in basis}
    # 1. exact symbol
    t = {}
    for nm, a, b in zip(names, ca, cb):
        t = add(t, scale(symmap[nm], Q(a, b)))
    n5 = json.loads((HERE/'cobracket_e.json').read_text())[name.upper()]
    ts = WS.symbol(WS.e_expression(which))
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        ts = add(ts, scale(g_symbol(ZERO, li32_word(x, y), ONE), -c))
    assert t == ts, 'symbol mismatch'
    print(f"1. EXACT symbol of R({name}) equals symbol({name.upper()}) - symbol(N): PASS", flush=True)
    # 2. cobracket
    for nm, a, b in zip(names, ca, cb):
        if recmap[nm][0] == 'li' and recmap[nm][1] == 5:
            z = recmap[nm][2]
            assert not cred(delta_word(ZERO, (ZERO, ZERO, ZERO, ZERO, ONE/z), ONE)).d, nm
    print("2a. every Li5 atom in R reduces to 0 in the (3,2) cobracket: PASS", flush=True)
    prod = P()
    import w5_symbol as W5
    for w in W5.shuffles((ZERO, ONE/B), (ZERO, ZERO, ONE/B)):
        prod += delta_word(ZERO, w, ONE)
    assert not cred(prod).d
    print('2b. shuffle-product cobracket control: PASS', flush=True)
    right = P()
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        right += delta_word(ZERO, li32_word(x, y), ONE) * c
    tgt = C.P()
    for (w, k), c in WS.e_expression(which).items():
        if k == 0:
            tgt += (delta_word(ONE, w, A) - delta_word(ONE, w, B)) * c
    delta = cred(right - tgt)
    assert not delta.d
    print(f"2c. EXACT (3,2) cobracket: reduce(N({name}) - target({name})) = 0: PASS", flush=True)
    # 3. value at 1400 digits
    ctx.prec = math.ceil(1600*math.log2(10)); mp.mp.dps = 1530
    total = mp.mpf(0)
    for nm, a, b in zip(names, ca, cb):
        v = evaluate(recmap[nm], 1500)
        total += (mp.mpf(a.numerator)/a.denominator + mp.sqrt(2)*mp.mpf(b.numerator)/b.denominator) \
            * mp.mpf(v.mid().str(1520, radius=False))
    n5v = mp.mpf(0)
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        v = evaluate(('mpl', 3, 2, x, y), 1500)
        n5v += qnum(c)*mp.mpf(v.mid().str(1520, radius=False))
    total += n5v
    ev = mp.mpf(json.loads((HERE/f'{name}_5800.json').read_text())['value_full'])
    res = total - ev
    print(f"3. R({name}) + N({name}) - {name}(5800 bits): {mp.nstr(res, 8)}  (expect ~1e-1500)", flush=True)
    print(f'ALL {name} VERIFICATION CHECKS DONE')
    with open(HERE/'progress.log','a') as f:
        f.write(f'step_K verify {name}{suffix}: symbol EXACT PASS; cobracket EXACT PASS; value residual {mp.nstr(res, 8)}\n')

main()
