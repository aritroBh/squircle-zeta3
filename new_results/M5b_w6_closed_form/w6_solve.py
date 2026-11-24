"""W6 - N6: exact symbol stage.
Atoms: weight-6 products (weight-5 factors = Li5 at R5_min3 args, z5, Lchi, 4 Li32 atoms)
+ Li6 at S-unit args. Target = symbol(W6) - symbol(N6). Saturated kernel by HNF.
Usage: -B w6_solve.py [extra_li6_index ...]  (pool extras beyond the base 10)
"""
import itertools, json, re, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from w4_symbol import Q, R, A, B, ban_float, self_check, g_symbol, add, scale, li_symbol
import w5_symbol as W
import w5pp_symbol as WS
from m5b_common import (ZERO, ONE, log, pool_zs, atoms_w5, li32_word, symbol_solve,
                        saturated_kernel, product_atoms, W4_SEVEN)

WORDS6 = list(itertools.product(range(2), repeat=6))

def parse_qs2(s):
    s = s.strip().strip('()').replace(' ', '')
    if 'sqrt2' not in s:
        return Q(F(s))
    s = s.replace('*sqrt2', 'S')
    p, q = F(0), F(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if t.endswith('S'):
            c = t[:-1]
            if c in ('', '+'):
                c = '1'
            elif c == '-':
                c = '-1'
            q += F(c)
        else:
            p += F(t)
    return Q(p, q)

def n6_terms():
    d = json.loads((HERE/'n6_solution.json').read_text())
    out = []
    for t in d['atoms']:
        x = parse_qs2(t['x']); y = parse_qs2(t['y'])
        c = Q(F(t['coeff_p']), F(t['coeff_q']))
        out.append((t['m'], t['n'], x, y, c))
    return out

def mpl_word(m, n, x, y):
    return (ZERO,)*(m-1) + (ONE/x,) + (ZERO,)*(n-1) + (ONE/(x*y),)

def target_symbol6():
    t = WS.symbol(WS.expression_from([(3, 2), (4, 1), (5, 0)]))
    for m, n, x, y, c in n6_terms():
        t = add(t, scale(g_symbol(ZERO, mpl_word(m, n, x, y), ONE), -c))
    return t

R5_MIN3 = [B, -B, -B*B, R/2, Q(F(1, 2)), ONE-R/2, 2-R,
           Q(F(1, 2))-R/2, Q(F(1, 2))-R/4, Q(F(1, 2))+R/4]
LI32_ATOMS = [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE)]

def factor_groups6():
    g5 = [(f'Li5({z})', W.li_symbol(5, z), ('li', 5, z)) for z in R5_MIN3]
    g5.append(('z5', {}, ('zeta', 5)))
    g5.append(('Lchi8_5', {}, ('lchi',)))
    for x, y in LI32_ATOMS:
        g5.append((f'Li32({x},{y})', g_symbol(ZERO, li32_word(x, y), ONE), ('mpl', 3, 2, x, y)))
    return {
        1: [('lam', W.log_symbol(Q(2)), ('log', Q(2))),
            ('L', W.log_symbol(A), ('log', A))],
        2: [('pi2', {}, ('pi2',)), ('Li2(b)', W.li_symbol(2, B), ('li', 2, B))],
        3: [('z3', {}, ('zeta', 3))] + [(f'Li3({z})', W.li_symbol(3, z), ('li', 3, z))
                                        for z in [B, -B, R/2]],
        4: [(f'Li4({z})', W.li_symbol(4, z), ('li', 4, z)) for z in W4_SEVEN],
        5: g5,
    }

BASE_LI6 = W4_SEVEN + [Q(F(1, 2))-R/2, Q(F(1, 2))-R/4, Q(F(1, 2))+R/4]  # + -b/2, (2-+sqrt2)/4

def atoms_w6(li6_args):
    out = product_atoms(6, factor_groups6())
    for z in li6_args:
        out.append((f'Li6({z})', li_symbol(6, z), ('li', 6, z)))
    return out

if __name__ == '__main__':
    ban_float(); self_check()
    extra_idx = [int(i) for i in sys.argv[1:]]
    pool = pool_zs()
    li6 = list(BASE_LI6)
    for i in extra_idx:
        if pool[i] not in li6:
            li6.append(pool[i])
    print('Li6 arguments:', [str(z) for z in li6], flush=True)
    basis = atoms_w6(li6)
    print('atom count:', len(basis), flush=True)
    t = target_symbol6()
    sol = symbol_solve(basis, t, WORDS6)
    if not sol['ok']:
        log('W6 symbol stage: target NOT in span')
        sys.exit(1)
    import math
    mx = max(max(abs(c.numerator), 1) for c in sol['xp']+sol['xq'])
    dn = math.lcm(*[c.denominator for c in sol['xp']+sol['xq']])
    print(f'rref particular: max|num|={mx} lcm-den={dn}; kernel dim={len(sol["kernel"])}', flush=True)
    K, rank = saturated_kernel(basis, WORDS6)
    names = [b[0] for b in basis]
    tag = 'w6' if not extra_idx else 'w6x' + '_'.join(map(str, extra_idx))
    data = dict(names=names, rank=rank, li6_args=[str(z) for z in li6],
                xp=[str(c) for c in sol['xp']], xq=[str(c) for c in sol['xq']],
                Ksat=K)
    (HERE/f'{tag}_symbol_solve.json').write_text(json.dumps(data)+'\n')
    log(f'W6 symbol stage done: rank={rank}, kernel_sat_dim={len(K)}, particular max|num|={mx}, lcm={dn}')
