"""WT5 step G: E_i - N(E_i): exact symbol stage. Target = symbol(E_i) - symbol(N_i) in the span of
products + Li5(16 pool S-units) + z5. Saturated kernel by HNF. Saves e{i}_symbol_solve.json.
Usage: -B wt5_G_solve.py <e1|e2|e3>"""
import itertools, json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w4_symbol import Q, R, A, B, ban_float, self_check, g_symbol, add, scale
import w5_symbol as W
import wt5_E_pairings as WS
from m5b_common import (ZERO, ONE, log, pool_zs, atoms_w5, li32_word, symbol_solve,
                        saturated_kernel)

WORDS5 = list(itertools.product(range(2), repeat=5))

def n_terms(name):
    d = json.loads((HERE/'cobracket_e.json').read_text())
    out = []
    for t in d[name.upper()]['correction']:
        x = Q(*map(F, t['x'])); y = Q(*map(F, t['y'])); c = Q(*map(F, t['coefficient']))
        out.append((x, y, c))
    return out

def target_symbol(name):
    which = {'e1': 'D1', 'e2': 'D2', 'e3': 'D3'}[name]
    t = WS.symbol(WS.e_expression(which))
    for x, y, c in n_terms(name):
        t = add(t, scale(g_symbol(ZERO, li32_word(x, y), ONE), -c))
    return t

if __name__ == '__main__':
    name = sys.argv[1]
    ban_float(); self_check()
    zs = pool_zs()
    basis = atoms_w5(zs)
    print('atom count:', len(basis), flush=True)
    t = target_symbol(name)
    # extra diagnostic for e3: support of the symbol by letter
    if name == 'e3':
        letters = sorted({x for w in t for x in w})
        only0 = all(set(w) <= {0} for w in t)
        print('E3 target symbol letter support:', letters, ' only letter 0 (log2):', only0, flush=True)
    sol = symbol_solve(basis, t, WORDS5)
    assert sol['ok'], 'target symbol not in span -- enlarge Li5 set'
    mx = max(max(abs(c.numerator), 1) for c in sol['xp']+sol['xq'])
    import math
    dn = math.lcm(*[c.denominator for c in sol['xp']+sol['xq']])
    print(f'rref particular: max|num|={mx} lcm-den={dn}; kernel dim={len(sol["kernel"])}', flush=True)
    K, rank = saturated_kernel(basis, WORDS5)
    names = [b[0] for b in basis]
    data = dict(names=names, rank=rank,
                xp=[str(c) for c in sol['xp']], xq=[str(c) for c in sol['xq']],
                kernel_rref=[[str(c) for c in v] for v in sol['kernel']],
                Ksat=K)
    (HERE/f'{name}_symbol_solve.json').write_text(json.dumps(data)+'\n')
    log(f"{name} symbol stage done: rank={rank}, kernel_sat_dim={len(K)}, particular max|num|={mx}, lcm={dn}")
