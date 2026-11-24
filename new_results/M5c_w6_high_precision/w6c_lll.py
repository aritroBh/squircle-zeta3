"""M5c W6-N6 kernel LLL at <digits> on 3500-digit values (179+ columns).
Usage: -B w6c_lll.py <digits> [tag]   tag='' -> w6_symbol_solve.json + w6_{leaves,atomvals}_3500.pkl + w6_11800.json
                                       tag='w6x14' -> w6x14_symbol_solve.json + w6x14_atomvals_3500.pkl etc.
Planted-relation control first. Writes w6c_lll_<tag>_<digits>.json. Prints top-8 rows by |residual|.
"""
import json, math, pickle, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
from w4_symbol import Q, R, A, B
from m5b_common import log, qnum
from w6_solve import n6_terms
from w6_values import leaf_recipes

DIG_VALS = 3500
W6_JSON = 'w6_11800.json'

def set_prec(digits):
    mp.mp.dps = digits + 120

def build(tag, w6_json=W6_JSON):
    set_prec(DIG_VALS)
    solf = 'w6_symbol_solve.json' if not tag else f'{tag}_symbol_solve.json'
    sol = json.loads((HERE/solf).read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    suf = '' if not tag else '_' + tag.replace('w6x','x')
    avf = HERE/f'w6_atomvals_{DIG_VALS}{suf}.pkl'
    lvf = HERE/f'w6_leaves_{DIG_VALS}{suf}.pkl'
    if not lvf.exists():
        lvf = HERE/f'w6_leaves_{DIG_VALS}.pkl'
    vals = [mp.mpf(v) for v in pickle.loads(avf.read_bytes())]
    recs = leaf_recipes()
    saved = pickle.loads(lvf.read_bytes())
    if isinstance(saved, dict):
        saved = list(saved.values())
    leaves = {r: mp.mpf(v) for r, v in zip(recs, saved)}
    n6v = mp.mpf(0)
    for m, n, x, y, c in n6_terms():
        n6v += qnum(c)*leaves[('mpl', m, n, x, y)]
    w6 = mp.mpf(json.loads((HERE/w6_json).read_text())['value_full'])
    pval = mp.fsum(qnum(Q(a, b_))*v for a, b_, v in zip(xp0, xq0, vals))
    target = w6 - n6v - pval
    kv = [mp.fsum(F(c)*v for c, v in zip(krow, vals)) for krow in K]
    log(f'build[{tag or "w6"}]: W6-N6-particular target = {mp.nstr(target, 30)}; kernel dim={len(K)}; atoms={len(names)}')
    return dict(names=names, K=K, kv=kv, target=target)

def col_values(d, kept):
    r2 = mp.sqrt(2)
    return [d['target']] + [(r2 if j % 2 else 1)*d['kv'][j//2] for j in kept]

def planted():
    set_prec(600)
    n = 3
    C = mp.mpf(10)**500
    plant = [mp.mpf(1), mp.sqrt(2), 1+mp.sqrt(2)]
    rows = [[int(i == j) for j in range(n)]+[int(mp.nint(C*v))] for i, v in enumerate(plant)]
    red = flint.fmpz_mat(rows).lll(gram='exact')
    ok = any([int(red[i, j]) for j in range(n)] in ([1, 1, -1], [-1, -1, 1]) for i in range(n))
    assert ok
    log('planted relation control: PASS')

def lll(vals, digits, label):
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i == j) for j in range(n)]+[int(mp.nint(C*v))] for i, v in enumerate(vals)]
    log(f'LLL start {label}: rows={n}, digits={digits}, gram=exact')
    t0 = time.monotonic()
    reduced = flint.fmpz_mat(rows).lll(gram='exact')
    dt = time.monotonic()-t0
    set_prec(digits+200)
    relations = []
    diag = []
    for i in range(n):
        vec = [int(reduced[i, j]) for j in range(n)]
        mx = max(map(abs, vec))
        residual = mp.fsum(mp.mpf(c)*v for c, v in zip(vec, vals))
        diag.append((mx, residual, i))
        if abs(residual) < mp.mpf(10)**(-digits-60):
            relations.append(dict(vec=vec, maxcoef=mx, residual=mp.nstr(residual, 12)))
    log(f'LLL done {label}: seconds={dt:.2f}; small-residual relations={len(relations)}')
    diag.sort(key=lambda t: abs(t[1]))
    for mx, residual, i in diag[:8]:
        log(f'  row {i}: maxcoef={mx} residual={mp.nstr(residual, 8)}')
    return relations, [(str(mx), mp.nstr(res, 8), i) for mx, res, i in diag[:8]]

def main():
    digits = int(sys.argv[1])
    tag = sys.argv[2] if len(sys.argv) > 2 else ''
    planted()
    d = build(tag)
    nk = len(d['kv'])
    kept = list(range(2*nk))
    set_prec(digits+120)
    vals = [mp.mpf(v) for v in col_values(d, kept)]
    relations, top = lll(vals, digits, f'W6c symbol-kernel enlarged basis[{tag or "w6"}] @{digits}')
    out = dict(digits=digits, tag=tag, kept=kept, relations=relations, top_rows=top)
    (HERE/f'w6c_lll_{tag}_{digits}.json'.replace("__","_")).write_text(json.dumps(out, indent=2)+'\n')
    nrel = [r for r in relations if r['vec'][0] != 0]
    log(f'W6c LLL[{tag or "w6"}] @{digits}: target-involving relations={len(nrel)}')
    for r in nrel:
        log(f'  c0={r["vec"][0]} maxcoef={r["maxcoef"]} residual={r["residual"]}')

if __name__ == '__main__':
    main()
