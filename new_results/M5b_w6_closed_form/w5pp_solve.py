"""W5'' - N5'': exact symbol stage. Target = symbol(W5'') - symbol(N5'') in the span of
products + Li5(16 pool S-units) + z5. Saturated kernel by HNF. Saves w5pp_symbol_solve.json."""
import itertools, json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w4_symbol import Q, R, A, B, ban_float, self_check, g_symbol, add, scale
import w5_symbol as W
import w5pp_symbol as WS
from m5b_common import (ZERO, ONE, log, pool_zs, atoms_w5, li32_word, symbol_solve,
                        saturated_kernel)

WORDS5 = list(itertools.product(range(2), repeat=5))

def n5pp_terms():
    d = json.loads((HERE/'cobracket5pp.json').read_text())
    out = []
    for t in d['correction']:
        x = Q(*map(F, t['x'])); y = Q(*map(F, t['y'])); c = Q(*map(F, t['coefficient']))
        out.append((x, y, c))
    return out

def target_symbol5pp():
    t = WS.symbol(WS.w5pp_expression())
    for x, y, c in n5pp_terms():
        t = add(t, scale(g_symbol(ZERO, li32_word(x, y), ONE), -c))
    return t

if __name__ == '__main__':
    ban_float(); self_check()
    zs = pool_zs()
    print('Li5 arguments (pool):', [str(z) for z in zs], flush=True)
    basis = atoms_w5(zs)
    print('atom count:', len(basis), flush=True)
    t = target_symbol5pp()
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
    (HERE/'w5pp_symbol_solve.json').write_text(json.dumps(data)+'\n')
    log(f"W5'' symbol stage done: rank={rank}, kernel_sat_dim={len(K)}, particular max|num|={mx}, lcm={dn}")
