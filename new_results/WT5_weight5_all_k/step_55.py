base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step H: E_i closed form: atom values (cached), kv kernel values, LLL with dependency dropping.
Usage: -B wt5_H_lll.py <e1|e2|e3> <digits>. Values cache: w5pp_vals_1450.pkl (same 71-atom basis)."""
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
import w5_symbol as W
from w5_search import evaluate
from m5b_common import ZERO, ONE, log, pool_zs, atoms_w5, li32_word, qnum, lll_relations

DIG_VALS = 1450

def set_prec(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120

def build(name, digits):
    set_prec(digits)
    sol = json.loads((HERE/f'{name}_symbol_solve.json').read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    zs = pool_zs()
    basis = atoms_w5(zs)
    assert [b[0] for b in basis] == names
    cache = HERE/f'w5pp_vals_{DIG_VALS}.pkl'
    if cache.exists():
        vals = pickle.loads(cache.read_bytes())
    else:
        set_prec(DIG_VALS)
        vals = []
        t0 = time.monotonic()
        for nm, _, recipe in basis:
            v = evaluate(recipe, DIG_VALS+40)
            vals.append(mp.mpf(v.mid().str(DIG_VALS+110, radius=False)))
        cache.write_bytes(pickle.dumps(vals))
        log(f'atom values computed ({time.monotonic()-t0:.1f}s)')
    set_prec(digits)
    vals = [mp.mpf(v) for v in vals]
    # N(E_i) value
    n5 = json.loads((HERE/'cobracket_e.json').read_text())[name.upper()]
    n5v = mp.mpf(0)
    for t in n5['correction']:
        x, y, c = (Q(*map(F, t[k])) for k in ['x', 'y', 'coefficient'])
        v = evaluate(('mpl', 3, 2, x, y), DIG_VALS+40)
        n5v += qnum(c)*mp.mpf(v.mid().str(DIG_VALS+110, radius=False))
    ev = mp.mpf(json.loads((HERE/f'{name}_5800.json').read_text())['value_full'])
    pval = mp.fsum(qnum(Q(a, b_))*v for a, b_, v in zip(xp0, xq0, vals))
    target = ev - n5v - pval
    kv = [mp.fsum(F(c)*v for c, v in zip(krow, vals)) for krow in K]
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    return dict(names=names, vals=vals, K=K, kv=kv, target=target, Lchi=Lchi,
                ev=ev, n5v=n5v, pval=pval)

def run(name, digits):
    d = build(name, digits)
    kv, Lchi, target = d['kv'], d['Lchi'], d['target']
    r2 = mp.sqrt(2)
    kept = list(range(2*len(kv)+2))   # kv0,sqrt2kv0,...,L,sqrt2 L
    def values():
        return [target]+[(r2 if j % 2 else 1)*(kv[j//2] if j//2 < len(kv) else Lchi) for j in kept]
    relations = lll_relations(values(), digits, f"{name} symbol-kernel enlarged basis")
    dependencies = []
    while any(r['vec'][0] == 0 for r in relations):
        deps = [r for r in relations if r['vec'][0] == 0]
        dependencies.append(dict(kept=kept.copy(), relations=deps))
        mat, rank = flint.fmpq_mat([list(reversed(r['vec'][1:])) for r in deps]).rref()
        drop = {len(kept)-1-next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)}
        kept = [j for i, j in enumerate(kept) if i not in drop]
        log(f'Removed {len(drop)} independent numerical basis dependencies; retained scalar atoms={len(kept)}')
        relations = lll_relations(values(), digits, f"{name} independent symbol-kernel enlarged basis")
    data = dict(name=name, digits=digits, kept=kept, relations=relations, dependencies=dependencies,
                target=mp.nstr(target, 50))
    (HERE/f'{name}_lll_{digits}.json').write_text(json.dumps(data, indent=2)+'\\n')
    print(json.dumps(dict(name=name, digits=digits, n_rel=len(relations),
                          rels=[(r['vec'], r['residual']) for r in relations]), indent=1), flush=True)

if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1100)
'''
open(W+'wt5_H_lll.py','w').write(script)
print('written')
# verify the atom value cache matches the basis
import pickle, json
W2=W
vals=pickle.load(open(W2+'w5pp_vals_1450.pkl','rb'))
print('cache len', len(vals), type(vals[0]))
