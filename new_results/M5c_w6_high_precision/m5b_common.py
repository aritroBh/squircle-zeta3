"""M5b shared: weight-5/6 atom builders, symbol solve, saturated kernel (HNF), values, LLL.
Atoms carry (name, symbol, recipe) with recipes for w5_search.evaluate.
"""
import itertools, json, math, time
from fractions import Fraction as F
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, acb, ctx
import mpmath as mp
from w4_symbol import Q, R, A, B
import w4_symbol as W4
import w5_symbol as W

ZERO, ONE = Q(0), Q(1)

def log(s):
    print(s, flush=True)
    import datetime
    with (HERE/'progress.log').open('a') as f:
        f.write(f'{datetime.datetime.now().isoformat()} {s}\n')

def pool_zs():
    d = json.loads((HERE/'pool.json').read_text())
    return [Q(F(e['p']), F(e['q'])) for e in d['pool']]

W4_SEVEN = [B, -B, -B*B, R/2, Q(F(1,2)), ONE-R/2, 2-R]  # W4's seven arguments

def factor_groups5():
    return {
        1: [('lam', W.log_symbol(Q(2)), ('log', Q(2))),
            ('L', W.log_symbol(A), ('log', A))],
        2: [('pi2', {}, ('pi2',)), ('Li2(b)', W.li_symbol(2, B), ('li', 2, B))],
        3: [('z3', {}, ('zeta', 3))] + [(f'Li3({z})', W.li_symbol(3, z), ('li', 3, z))
                                        for z in [B, -B, R/2]],
        4: [(f'Li4({z})', W.li_symbol(4, z), ('li', 4, z)) for z in W4_SEVEN],
    }

def product_atoms(weight, groups):
    out = []
    def partitions(n, mx=None):
        if n == 0:
            yield []
        for k in (mx or sorted(groups, reverse=True)):
            if k <= n:
                for rest in partitions(n-k, [x for x in groups if x <= k]):
                    yield [k]+rest
    seen = set()
    for part in partitions(weight):
        part = tuple(sorted(part, reverse=True))
        if part in seen:
            continue
        seen.add(part)
        choices = [list(itertools.combinations_with_replacement(groups[k], part.count(k)))
                   for k in sorted(set(part))]
        for blocks in itertools.product(*choices):
            atoms = sum(blocks, ())
            s = {(): ONE}
            for _, v, _ in atoms:
                s = W.shuffle(s, v)
            out.append(('*'.join(a[0] for a in atoms), s, ('product', tuple(a[2] for a in atoms))))
    return out

def li32_word(x, y):
    return (ZERO, ZERO, ONE/x, ZERO, ONE/(x*y))

def depth2_w5_atoms():
    return [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE)]

def atoms_w5(li5_args):
    groups = factor_groups5()
    out = product_atoms(5, groups)
    for z in li5_args:
        out.append((f'Li5({z})', W.li_symbol(5, z), ('li', 5, z)))
    out.append(('z5', {}, ('zeta', 5)))
    return out

# ---------------- exact symbol solve over Q(sqrt2): atoms rational, target p+sqrt2 q
def symbol_solve(basis, target, words):
    assert all(c.q == 0 for _, s, _ in basis for c in s.values())
    n = len(basis)
    rows = [[flint.fmpq(str(s.get(w, ZERO).p)) for _, s, _ in basis]
            + [flint.fmpq(str(target.get(w, ZERO).p)), flint.fmpq(str(target.get(w, ZERO).q))]
            for w in words]
    t0 = time.monotonic()
    mat, rank = flint.fmpq_mat(rows).rref()
    pivots = [next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)]
    ok = all(j < n for j in pivots)
    log(f'symbol solve: atoms={n} words={len(words)} rank={rank} target-in-span={ok} ({time.monotonic()-t0:.1f}s)')
    if not ok:
        return dict(ok=False, rank=rank)
    free = [j for j in range(n) if j not in pivots]
    xp = [F(0)]*n; xq = [F(0)]*n
    for i, p in enumerate(pivots):
        xp[p] = F(str(mat[i, n])); xq[p] = F(str(mat[i, n+1]))
    kernel = []
    for j in free:
        v = [F(0)]*n
        v[j] = F(1)
        for i, p in enumerate(pivots):
            v[p] = -F(str(mat[i, j]))
        kernel.append(v)
    # independent reconstruction check
    s = dict(target)
    for j in range(n):
        c = Q(xp[j], xq[j])
        if c:
            s = W.add(s, W.scale(basis[j][1], -c))
    assert not s, 'particular reconstruction failed'
    return dict(ok=True, rank=rank, xp=xp, xq=xq, kernel=kernel, pivots=pivots, free=free)

# ---------------- saturated integer kernel via HNF of transpose
def saturated_kernel(basis, words):
    from math import gcd, lcm
    n = len(basis)
    cols = []
    for name, s, _ in basis:
        assert all(c.q == 0 for c in s.values()), name
        cols.append([s.get(w, ZERO).p for w in words])
    D2 = 1
    for col in cols:
        for c in col:
            D2 = lcm(D2, c.denominator)
    M2 = [[int(c*D2) for c in col] for col in cols]   # atoms x words
    rows = [[M2[ai][wi] for ai in range(n)] for wi in range(len(words))]  # words x atoms
    Amat = flint.fmpz_mat(rows)
    t0 = time.monotonic()
    H, U = Amat.transpose().hnf(transform=True)
    rank = sum(1 for i in range(H.nrows()) if any(H[i, j] for j in range(H.ncols())))
    Ksat = [[int(U[i, j]) for j in range(n)] for i in range(rank, n)]
    for krow in Ksat:
        for wi in range(len(words)):
            assert sum(F(M2[ai][wi])*F(krow[ai]) for ai in range(n)) == 0
    snf = flint.fmpz_mat(Ksat).snf()
    diag = [int(snf[i, i]) for i in range(min(snf.nrows(), snf.ncols()))]
    log(f'HNF saturation: {time.monotonic()-t0:.1f}s rank={rank} kernel_dim={len(Ksat)} snf_invariants={sorted(set(diag))}')
    t0 = time.monotonic()
    red = flint.fmpz_mat(Ksat).lll(gram='exact')
    K = [[int(red[i, j]) for j in range(red.ncols())] for i in range(red.nrows())]
    log(f'LLL-reduce saturated kernel: {time.monotonic()-t0:.1f}s; max row entry {max(max(abs(i) for i in r) for r in K)}')
    return K, rank

# ---------------- values
def qnum(c):
    return mp.mpf(c.p.numerator)/c.p.denominator + mp.sqrt(2)*mp.mpf(c.q.numerator)/c.q.denominator

def atom_values(basis, digits):
    from w5_search import evaluate
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120
    out = []
    t0 = time.monotonic()
    for name, _, recipe in basis:
        v = evaluate(recipe, digits+40)
        out.append(mp.mpf(v.mid().str(digits+110, radius=False)))
        log(f'  value {name} ({time.monotonic()-t0:.1f}s)')
    return out

def lll_relations(vals, digits, label):
    import signal
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i == j) for j in range(n)]+[int(mp.nint(C*v))] for i, v in enumerate(vals)]
    log(f'LLL start {label}: rows={n}, digits={digits}, gram=exact')
    signal.alarm(3300)
    t0 = time.monotonic()
    reduced = flint.fmpz_mat(rows).lll(gram='exact')
    signal.alarm(0)
    relations = []
    for i in range(n):
        vec = [int(reduced[i, j]) for j in range(n)]
        mx = max(map(abs, vec))
        residual = mp.fsum(c*v for c, v in zip(vec, vals))
        if abs(residual) < mp.mpf(10)**(-digits-30) and mx < mp.mpf(10)**(mp.mpf(digits)/n):
            relations.append(dict(vec=vec, maxcoef=mx, residual=mp.nstr(residual, 12)))
    log(f'LLL done {label}: seconds={time.monotonic()-t0:.2f}; small relations={len(relations)}')
    return relations
