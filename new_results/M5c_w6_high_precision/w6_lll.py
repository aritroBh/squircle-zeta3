"""W6-N6 kernel LLL: planted test, dependency hunt, final relation.
Usage: -B w6_lll.py <digits> [depdigits]
Columns: [target] + kv_j, sqrt2 kv_j for the 89 saturated-kernel combos.
"""
import json, math, pickle, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, acb, ctx
import mpmath as mp
from w4_symbol import Q, R, A, B
from m5b_common import log, qnum, lll_relations
from w6_solve import n6_terms, BASE_LI6
from w6_values import leaf_recipes

DIG_VALS = 2400

def set_prec(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120

def build(digits):
    set_prec(digits)
    sol = json.loads((HERE/'w6_symbol_solve.json').read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    vals = [mp.mpf(v) for v in pickle.loads((HERE/f'w6_atomvals_{DIG_VALS}.pkl').read_bytes())]
    recs = leaf_recipes()
    saved = pickle.loads((HERE/f'w6_leaves_{DIG_VALS}.pkl').read_bytes())
    if isinstance(saved, dict):
        saved = list(saved.values())
    leaves = {r: mp.mpf(v) for r, v in zip(recs, saved)}
    n6v = mp.mpf(0)
    for m, n, x, y, c in n6_terms():
        n6v += qnum(c)*leaves[('mpl', m, n, x, y)]
    w6 = mp.mpf(json.loads((HERE/'w6_8000.json').read_text())['value_full'])
    pval = mp.fsum(qnum(Q(a, b_))*v for a, b_, v in zip(xp0, xq0, vals))
    target = w6 - n6v - pval
    kv = [mp.fsum(F(c)*v for c, v in zip(krow, vals)) for krow in K]
    log(f'build: W6-N6-particular target = {mp.nstr(target, 30)}')
    return dict(names=names, K=K, kv=kv, target=target, w6=w6, n6v=n6v, pval=pval)

def col_values(d, kept):
    r2 = mp.sqrt(2)
    return [d['target']] + [(r2 if j % 2 else 1)*d['kv'][j//2] for j in kept]

def main():
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
    depdigits = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    d = build(digits)
    nk = len(d['kv'])
    # planted relation control
    set_prec(600)
    plant = [mp.mpf(1), mp.sqrt(2), 1+mp.sqrt(2)]
    planted = lll_relations(plant, 500, 'planted relation')
    assert any(r['vec'] in [[1, 1, -1], [-1, -1, 1]] for r in planted), planted
    log('planted relation control: PASS')
    # dependency hunt on kernel columns only
    set_prec(depdigits+100)
    allcols = list(range(2*nk))
    depvals = [mp.mpf(v) for v in col_values(d, allcols)[1:]]
    deps = lll_relations(depvals, depdigits, f'W6 dependency hunt (kernel columns only)')
    # verify candidates at full precision
    set_prec(DIG_VALS)
    fullvals = [mp.mpf(v) for v in col_values(d, allcols)[1:]]
    exact = []
    for r in deps:
        res = mp.fsum(F(c)*v for c, v in zip(r['vec'], fullvals))
        if abs(res) < mp.mpf(10)**(-(DIG_VALS-60)):
            exact.append((r['vec'], mp.nstr(res, 6)))
    log(f'dependency candidates={len(deps)}; verified exact at {DIG_VALS} digits: {len(exact)}')
    for v, res in exact:
        log(f'  dep maxcoef={max(map(abs,v))} residual1650={res}')
    kept = allcols
    if exact:
        mat, rank = flint.fmpq_mat([list(reversed(v)) for v, _ in exact]).rref()
        drop = {len(kept)-1-next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)}
        kept = [j for i, j in enumerate(kept) if i not in drop]
        log(f'dropped {len(drop)} dependent columns; retained {len(kept)}')
    # final LLL
    set_prec(digits+100)
    relations = lll_relations(col_values(d, kept), digits, 'W6 symbol-kernel enlarged basis')
    out = dict(digits=digits, depdigits=depdigits, kept=kept,
               exact_deps=[dict(vec=v, residual=r) for v, r in exact],
               relations=relations, target=mp.nstr(d['target'], 60))
    (HERE/f'w6_lll_{digits}.json').write_text(json.dumps(out, indent=2)+'\n')
    nrel = [r for r in relations if r['vec'][0] != 0]
    log(f'W6 LLL @{digits}: target-involving relations={len(nrel)}')
    for r in nrel:
        log(f'  relation c0={r["vec"][0]} maxcoef={r["maxcoef"]} residual={r["residual"]}')
    print(json.dumps(dict(digits=digits, n_target_rel=len(nrel),
                          c0=[r['vec'][0] for r in nrel],
                          residual=[r['residual'] for r in nrel]), indent=1), flush=True)

if __name__ == '__main__':
    main()
