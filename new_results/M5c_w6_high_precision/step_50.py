extra='''"""M5c: leaf+atom values for extra Li6 pool args at 3560 digits.
Usage: -B w6c_values_extra.py <pool_index> [...]
Reads w6x<idx>_symbol_solve.json (must exist; run w6_solve.py first) for the atom basis,
base cache w6_leaves_3500.pkl, computes new ('li',6,z) leaves, writes
w6_leaves_extra_w6x<idxs>.pkl (dict recipe->mpf str) and w6_atomvals_3500_x<idxs>.pkl (list).
"""
import json, math, pickle, sys, time
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp
mp.mp.dps = 3700
from flint import arb, acb, ctx
from w4_symbol import Q, R, A, B
from m5b_common import log, pool_zs
from w6_solve import atoms_w6, BASE_LI6
from w6_values import leaf_recipes, numeric

ctx.prec = 12700

def li6_value(z, digits):
    v = acb(numeric(z)).polylog(6)
    assert v.imag.contains(0)
    return v.real

def main():
    idxs = [int(i) for i in sys.argv[1:]]
    pool = pool_zs()
    li6 = list(BASE_LI6)
    new_args = []
    for i in idxs:
        if pool[i] not in li6:
            li6.append(pool[i]); new_args.append(pool[i])
    tag = 'w6x' + '_'.join(map(str, idxs))
    saved = pickle.loads((HERE/'w6_leaves_3500.pkl').read_bytes())
    recs = leaf_recipes()
    leaves = {r: mp.mpf(v) for r, v in zip(recs, saved)}
    t0 = time.monotonic()
    extra = {}
    for z in new_args:
        rec = ('li', 6, z)
        v = li6_value(z, 3560)
        sv = v.mid().str(3640, radius=False)
        leaves[rec] = mp.mpf(sv)
        extra[rec] = sv
        log(f'  extra leaf li6({z}) rad={v.rad().str(6)} ({time.monotonic()-t0:.1f}s)')
    (HERE/f'w6_leaves_extra_{tag}.pkl').write_bytes(pickle.dumps(extra))
    basis = atoms_w6(li6)
    vals = []
    for name, _, recipe in basis:
        if recipe[0] == 'product':
            v = mp.mpf(1)
            for leaf in recipe[1]:
                v *= leaves[leaf]
            vals.append(v)
        else:
            vals.append(leaves[recipe])
    (HERE/f'w6_atomvals_3500_x{"_".join(map(str,idxs))}.pkl').write_bytes(
        pickle.dumps([mp.nstr(v, 3620) for v in vals]))
    log(f'extra atom values saved for {tag}: {len(vals)} atoms, new Li6 args={[str(z) for z in new_args]}')

if __name__ == '__main__':
    main()
'''
open('w6c_values_extra.py','w').write(extra)
import py_compile; py_compile.compile('w6c_values_extra.py', doraise=True)
print('w6c_values_extra.py ok')
