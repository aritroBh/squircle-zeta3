"""M4e: kernel-value LLL for R5_nice WITH L(chi8,5) as symbol-zero columns.

Atom set: M4d nice_symbol.json (S0; 65 atoms; rank 30; 35 rational kernel vectors).
Columns: [target] + kv_j, sqrt2*kv_j (j=0..34) + [Lchi, sqrt2*Lchi].
target = W5 - N5 - particular-value.  Genuine relation has vec[0] != 0.
Dependencies (vec[0]=0) are dropped iteratively (never dropping L or target).
"""
import json, math, pickle, signal, sys, time
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path

import flint
from flint import arb, ctx
import mpmath as mp

import w5_symbol as W
from w5_search import evaluate, mpl, numeric

HERE = Path(__file__).resolve().parent
ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B

def log(s):
    print(s, flush=True)

def qnum(c):
    return mp.mpf(c.p.numerator)/c.p.denominator + mp.sqrt(2)*mp.mpf(c.q.numerator)/c.q.denominator

def prim(vec):
    D = lcm(*[c.denominator for c in vec])
    ints = [int(c*D) for c in vec]
    g = 0
    for i in ints: g = gcd(g, abs(i))
    return [i//g for i in ints]

def lll_rows(vals, digits):
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i==j) for j in range(n)] + [int(mp.nint(C*v))] for i, v in enumerate(vals)]
    signal.alarm(2700)
    t0 = time.monotonic()
    red = flint.fmpz_mat(rows).lll(gram='exact')
    signal.alarm(0)
    out = []
    for i in range(n):
        vec = [int(red[i, j]) for j in range(n)]
        residual = mp.fsum(c*v for c, v in zip(vec, vals))
        out.append(dict(vec=vec, maxcoef=max(map(abs, vec)), residual=residual,
                        last=int(red[i, n])))
    log(f'LLL rows={n} digits={digits}: {time.monotonic()-t0:.1f}s')
    return out

def nice_score(colname):
    # lower = nicer = keep
    if colname in ('Lchi', 's2*Lchi'): return -10   # vital
    if colname.startswith('kv'):
        return 5
    return 9

def main(digits):
    ctx.prec = math.ceil((digits + 140) * math.log2(10))
    mp.mp.dps = digits + 120
    # planted relation test
    plant = lll_rows([mp.mpf(1), mp.sqrt(2), 1 + mp.sqrt(2)], digits)
    assert any(r['vec'] in [[1, 1, -1], [-1, -1, 1]] for r in plant)
    log('planted relation OK')
    # MPL control
    test = mpl(1, 1, B, ONE, digits + 20) - numeric(ONE - B).log()**2/2
    assert test.contains(0)
    log('MPL control OK')

    sym = json.loads((HERE/'nice_symbol.json').read_text())
    names = sym['names']
    xp = [F(s) for s in sym['xp']]
    xq = [F(s) for s in sym['xq']]
    ks = json.loads((HERE/'nice_ksat_s0.json').read_text())
    assert ks['names'] == names
    K = ks['K']
    log(f'kernel vectors={len(K)} max|entry|={max(max(map(abs,r)) for r in K)}')
    recipes = pickle.loads((HERE/'nice_atoms.pkl').read_bytes())
    assert [n for n, _ in recipes] == names

    t0 = time.monotonic()
    values = [evaluate(recipe, digits + 60) for _, recipe in recipes]
    maxrad = max(x.rad() for x in values)
    assert maxrad < arb(10)**(-digits - 30), maxrad
    vals = [mp.mpf(v.mid().str(digits + 110, radius=False)) for v in values]
    log(f'atom values: {time.monotonic()-t0:.1f}s maxrad={maxrad.str(4)}')
    r2 = mp.sqrt(2)

    cert = json.loads((HERE/'cobracket.json').read_text())
    fixedv = mp.mpf(0)
    for term in cert['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        v = evaluate(('mpl', 3, 2, x, y), digits + 60)
        fixedv += qnum(c)*mp.mpf(v.mid().str(digits + 110, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    assert digits + 80 < len(saved['value'])
    pval = mp.fsum((mp.mpf(x.numerator)/x.denominator + r2*mp.mpf(q.numerator)/q.denominator)*v
                   for x, q, v in zip(xp, xq, vals))
    target = mp.mpf(saved['value']) - fixedv - pval
    log(f'target (W5-N5-particular) ~ {mp.nstr(target, 12)}')

    kv = [mp.fsum(F(c)*v for c, v in zip(krow, vals)) for krow in K]
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])

    colnames = ['target'] + [f'kv{j}' for j in range(len(K))] + \
               [f's2kv{j}' for j in range(len(K))] + ['Lchi', 's2*Lchi']
    colvals = [target] + kv + [r2*v for v in kv] + [Lchi, r2*Lchi]
    kept = list(range(len(colvals)))
    tol = mp.mpf(10)**(-digits - 20)
    history = []
    for round_ in range(12):
        vals_ = [colvals[j] for j in kept]
        rows = lll_rows(vals_, digits)
        genuine = [r for r in rows if r['vec'][0] != 0 and abs(r['residual']) < tol]
        deps = [r for r in rows if r['vec'][0] == 0 and abs(r['residual']) < tol]
        genuine.sort(key=lambda r: r['maxcoef'])
        log(f'round {round_}: cols={len(kept)} genuine={len(genuine)} deps={len(deps)}'
            + (f' best maxcoef={genuine[0]["maxcoef"]} vec0={genuine[0]["vec"][0]}' if genuine else ''))
        history.append(dict(round=round_, kept=list(kept),
                            genuine=[dict(vec=r['vec'], maxcoef=r['maxcoef'],
                                          residual=mp.nstr(r['residual'], 8)) for r in genuine[:6]],
                            deps=[dict(vec=r['vec'], maxcoef=r['maxcoef'],
                                       residual=mp.nstr(r['residual'], 8)) for r in deps[:20]]))
        if genuine or not deps:
            (HERE/f'lchi_sat_search_{digits}.json').write_text(json.dumps(
                dict(digits=digits, colnames=colnames, K=K, history=history), indent=1)+'\n')
            return history
        # drop one column per dep: highest nice_score in support, never target(0)/L(last two)
        drop = set()
        for r in deps:
            sup = [i for i, c in enumerate(r['vec']) if c and i not in (0, len(r['vec'])-2, len(r['vec'])-1)]
            sup = [i for i in sup if kept[i] not in drop]
            if sup:
                worst = max(sup, key=lambda i: nice_score(colnames[kept[i]]))
                drop.add(kept[worst])
        log(f'dropping {len(drop)}: {[colnames[j] for j in sorted(drop)]}')
        kept = [j for j in kept if j not in drop]
    (HERE/f'lchi_sat_search_{digits}.json').write_text(json.dumps(
        dict(digits=digits, colnames=colnames, K=K, history=history), indent=1)+'\n')
    return history

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1200)
