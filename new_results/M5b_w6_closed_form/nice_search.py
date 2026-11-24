"""M4d: readable R5_nice. Adapted from M4c's w5_search.py (two-stage: exact symbol, then small LLL).

Differences from M4c:
 - weight-4 factors are W4's seven arguments {b,-b,-b^2,1/sqrt2,1/2,1-1/sqrt2,2-sqrt2};
 - Li5 arguments come from a NICE set S (start S0), not the 17 S-units;
 - the symbol stage reports membership and, on failure, the missing word directions;
 - the exact symbol solve keeps the full rational kernel, for coset nicification.
"""
import itertools
import json
import math
import pickle
import sys
from fractions import Fraction as F
from pathlib import Path

import flint
import w5_symbol as W

HERE = Path(__file__).resolve().parent
ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
WORDS = list(itertools.product(range(2), repeat=5))


def log(s):
    print(s, flush=True)


def factor_groups():
    return {
        1: [('lam', W.log_symbol(Q(2)), ('log', Q(2))),
            ('L', W.log_symbol(A), ('log', A))],
        2: [('pi2', {}, ('pi2',)), ('Li2(b)', W.li_symbol(2, B), ('li', 2, B))],
        3: [('z3', {}, ('zeta', 3))] + [(f'Li3({z})', W.li_symbol(3, z), ('li', 3, z))
                                        for z in [B, -B, R / 2]],
        4: [(f'Li4({z})', W.li_symbol(4, z), ('li', 4, z)) for z in
            [B, -B, -B * B, R / 2, Q(F(1, 2)), ONE - R / 2, 2 - R]],
    }


def products():
    groups = factor_groups()
    out = []
    for partition in [(1, 1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1), (3, 1, 1), (3, 2), (4, 1)]:
        choices = [list(itertools.combinations_with_replacement(groups[k], partition.count(k)))
                   for k in sorted(set(partition))]
        for blocks in itertools.product(*choices):
            atoms = sum(blocks, ())
            s = {(): ONE}
            for _, v, _ in atoms:
                s = W.shuffle(s, v)
            out.append(('*'.join(a[0] for a in atoms), s, ('product', tuple(a[2] for a in atoms))))
    return out


S0 = [B, -B, -B * B, R / 2, -R / 2, Q(F(1, 2)), ONE - R / 2, 2 - R, -ONE, B * B]
# M4c's remaining candidates, simplest first.
EXTRA = [-2 + 2 * R, B**4, Q(F(1, 2)) - R / 4, Q(F(1, 2)) + R / 4,
         Q(F(1, 2)) - R / 2, Q(F(1, 2)) - 3 * R / 8, -16 + 12 * R]


def atoms(li5_args):
    out = products()
    for z in li5_args:
        out.append((f'Li5({z})', W.li_symbol(5, z), ('li', 5, z)))
    out.append(('z5', {}, ('zeta', 5)))
    return out


def n5_correction():
    """The exact N5 from the freshly recomputed cobracket.json (must match M4c's)."""
    cert = json.loads((HERE / 'cobracket.json').read_text())
    fixed = []
    for term in cert['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x', 'y', 'coefficient'])
        word = (ZERO, ZERO, ONE / x, ZERO, ONE / (x * y))
        fixed.append((c, x, y, word))
    return fixed


def target_symbol():
    t = W.symbol()
    for c, x, y, word in n5_correction():
        t = W.add(t, W.scale(W.g_symbol(ZERO, word, ONE), -c))
    return t


def symbol_columns(basis):
    cols = []
    for name, s, recipe in basis:
        assert all(c.q == 0 for c in s.values()), name
        cols.append([F(str(s.get(w, ZERO).p)) for w in WORDS])
    return cols


def target_columns(t):
    return [F(str(t.get(w, ZERO).p)) for w in WORDS], [F(str(t.get(w, ZERO).q)) for w in WORDS]


def solve(li5_args, label=''):
    basis = atoms(li5_args)
    t = target_symbol()
    tp, tq = target_columns(t)
    cols = symbol_columns(basis)
    rows = [[str(c) for c in r] + [str(tp[i]), str(tq[i])] for i, r in enumerate(zip(*cols))]
    mat, rank = flint.fmpq_mat(rows).rref()
    n = len(basis)
    pivots = [next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)]
    ok = all(j < n for j in pivots)
    log(f'[{label}] atoms={n} symbol-rank={rank} target-in-span={ok}')
    if not ok:
        bad = [WORDS[i] for i in range(len(WORDS))
               if (tp[i] or tq[i]) and all(mat[r, n] == 0 and mat[r, n + 1] == 0 for r in range(rank))]
        miss = [i for i in range(len(WORDS)) if (tp[i] or tq[i])]
        # report target words that are zero in every basis column
        dead = [WORDS[i] for i in miss if all(c[i] == 0 for c in cols)]
        log(f'[{label}] target words zero in ALL atom symbols: {dead}')
        return dict(ok=False, rank=rank, basis=basis)
    free = [j for j in range(n) if j not in pivots]
    xp = [F(0)] * n
    xq = [F(0)] * n
    for i, p in enumerate(pivots):
        xp[p] = F(str(mat[i, n]))
        xq[p] = F(str(mat[i, n + 1]))
    kernel = []
    for j in free:
        v = [F(0)] * n
        v[j] = F(1)
        for i, p in enumerate(pivots):
            v[p] = -F(str(mat[i, j]))
        kernel.append(v)
    mx = max((max(abs(c.numerator), 1) for c in xp + xq), default=0)
    dn = math.lcm(*[c.denominator for c in xp + xq]) if xp else 1
    log(f'[{label}] free={len(free)} rref-particular max|num|={mx} lcm-den={dn}')
    return dict(ok=True, rank=rank, basis=basis, xp=xp, xq=xq, kernel=kernel,
                pivots=pivots, free=free, names=[b[0] for b in basis])


if __name__ == '__main__':
    W.ban_float()
    W.self_check()
    out = solve(S0, 'S0')
    if not out['ok']:
        for k, z in enumerate(EXTRA):
            out = solve(S0 + EXTRA[:k + 1], f'S0+extra[:{k+1}]={EXTRA[:k+1]}')
            if out['ok']:
                break
    if out['ok']:
        data = dict(names=out['names'], xp=[str(c) for c in out['xp']],
                    xq=[str(c) for c in out['xq']],
                    kernel=[[str(c) for c in v] for v in out['kernel']],
                    pivots=out['pivots'], free=out['free'], rank=out['rank'])
        (HERE / 'nice_symbol.json').write_text(json.dumps(data, indent=1) + '\n')
        recipes = [(b[0], b[2]) for b in out['basis']]
        (HERE / 'nice_atoms.pkl').write_bytes(pickle.dumps(recipes))
        log('saved nice_symbol.json and nice_atoms.pkl')
