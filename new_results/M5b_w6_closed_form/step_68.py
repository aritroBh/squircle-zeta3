import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""Verify W5'' = R5'' + N5'':
1. EXACT full point symbol of R5'' equals target (symbol(W5'')-symbol(N5'')).
2. EXACT (3,2) cobracket: each Li5 -> 0; reduce(N5'' - target(W5'')) = 0 (rerun M4c-style).
3. Value at 1350 digits vs W5pp_4800 (independent quadrature).
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
import w5pp_symbol as WS
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word
from w4_right_coproduct import P, PE
from w5_search import evaluate
from m5b_common import ZERO, ONE, pool_zs, atoms_w5, li32_word, qnum

def main():
    ban_float()
    d = json.loads((HERE/'w5pp_closed_form.json').read_text())
    names = d['names']; ca = [F(c) for c in d['ca']]; cb = [F(c) for c in d['cb']]
    aL, bL = F(d['aL']), F(d['bL'])
    zs = pool_zs()
    basis = atoms_w5(zs)
    symmap = {nm: s for nm, s, _ in basis}
    recmap = {nm: rec for nm, _, rec in basis}
    # 1. exact symbol
    t = {}
    for nm, a, b in zip(names, ca, cb):
        t = add(t, scale(symmap[nm], Q(a, b)))
    # target = symbol(W5'') - symbol(N5'')
    n5 = json.loads((HERE/'cobracket5pp.json').read_text())
    ts = WS.symbol(WS.w5pp_expression())
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        ts = add(ts, scale(g_symbol(ZERO, li32_word(x, y), ONE), -c))
    assert t == ts, 'symbol mismatch'
    print("1. EXACT symbol of R5'' equals symbol(W5'') - symbol(N5''): PASS", flush=True)
    # 2. cobracket
    for nm, a, b in zip(names, ca, cb):
        if recmap[nm][0] == 'li' and recmap[nm][1] == 5:
            z = recmap[nm][2]
            assert not cred(delta_word(ZERO, (ZERO, ZERO, ZERO, ZERO, ONE/z), ONE)).d, nm
    print("2a. every Li5 atom in R5'' reduces to 0 in the (3,2) cobracket: PASS", flush=True)
    # products reduce to 0 by construction of delta_word (kills shuffle products): spot check one
    prod = P()
    import w5_symbol as W5
    for w in W5.shuffles((ZERO, ONE/B), (ZERO, ZERO, ZERO, ONE/B)):
        prod += delta_word(ZERO, w, ONE)
    assert not cred(prod).d
    print('2b. shuffle-product cobracket control: PASS', flush=True)
    # N5'' vs target cobracket
    right = P()
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        right += delta_word(ZERO, li32_word(x, y), ONE) * c
    tgt = C.P()
    for (w, k), c in WS.w5pp_expression().items():
        if k == 0:
            tgt += (delta_word(ONE, w, A) - delta_word(ONE, w, B)) * c
    delta = cred(right - tgt)
    assert not delta.d
    print("2c. EXACT (3,2) cobracket: reduce(N5'' - target(W5'')) = 0: PASS", flush=True)
    # 3. value at 1350 digits
    ctx.prec = math.ceil(1500*math.log2(10)); mp.mp.dps = 1430
    total = mp.mpf(0)
    for nm, a, b in zip(names, ca, cb):
        v = evaluate(recmap[nm], 1400)
        total += (mp.mpf(a.numerator)/a.denominator + mp.sqrt(2)*mp.mpf(b.numerator)/b.denominator) \\
            * mp.mpf(v.mid().str(1420, radius=False))
    n5v = mp.mpf(0)
    for term in n5['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        v = evaluate(('mpl', 3, 2, x, y), 1400)
        n5v += qnum(c)*mp.mpf(v.mid().str(1420, radius=False))
    total += n5v
    w5pp = mp.mpf(json.loads((HERE/'w5pp_4800.json').read_text())['value_full'])
    res = total - w5pp
    print(f\"3. R5'' + N5'' - W5''(4800 bits): {mp.nstr(res, 8)}  (expect ~1e-1400)\", flush=True)
    print('ALL W5pp VERIFICATION CHECKS DONE')

main()
'''
open(os.path.join(m5b,'w5pp_verify.py'),'w').write(script)
import subprocess
PY='<project>/research/08_theory/M5b_w6_closed_form'
r=subprocess.run(['nice','-n','10',PY,'-B','w5pp_verify.py'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-2500:]); print(r.stderr[-900:])
