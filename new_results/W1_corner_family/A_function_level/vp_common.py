"""vp_common.py -- shared helpers for the vanishing-pair (VP) symbol-level systems.
Data format: output of a10b/a10c (terms N(X,Y), nice/bad Li4 columns, targets U,V symbols, ROOTS).
Functions are encoded as (s, k, e): f = s 2^k prod_r (c - r)^{e_r}, r in ROOTS.
"""
from fractions import Fraction as Fr
from collections import defaultdict
import flint
import mpmath as mp
from sym2 import delta22

mp.mp.dps = 30
GRID = None


def fval(f, x, ROOTS):
    s, k, e = f
    v = mp.mpf(s) * mp.mpf(2) ** k
    for r, ex in zip(ROOTS, e):
        v *= (x - mp.mpf(r.numerator) / r.denominator) ** ex
    return v


def path_ok(f, ROOTS, allow_above=False):
    global GRID
    if GRID is None:
        GRID = [mp.mpf(i) / 1000 for i in range(1, 708)] + [1 / mp.sqrt(2)]
    s, k, e = f
    for r, ex in zip(ROOTS, e):
        if ex < 0 and 0 < r and r * r <= Fr(1, 2):
            return False
    vals = [fval(f, x, ROOTS) for x in GRID]
    below = all(v < 1 - mp.mpf(10) ** -12 for v in vals)
    above = all(v > 1 + mp.mpf(10) ** -12 for v in vals)
    return below or (allow_above and above)


def at0(f, ROOTS):
    s, k, e = f
    if e[0] > 0:
        return Fr(0)
    if e[0] < 0:
        return 'inf'
    v = Fr(s) * Fr(2) ** k
    for r, ex in zip(ROOTS[1:], e[1:]):
        v *= (-r) ** ex
    return v


def cx(f):
    return sum(abs(x) for x in f[2]) + abs(f[1])


def term_ok(t, ROOTS, mode='lev2'):
    X, Y = t['X'], t['Y']
    if not (path_ok(X, ROOTS) and path_ok(Y, ROOTS)):
        return False
    a, b = at0(X, ROOTS), at0(Y, ROOTS)
    if mode == 'path':
        return True
    if a == 'inf' or b == 'inf':
        return False
    S = {a, b}
    if mode == 'lev2':
        return S <= {Fr(0), Fr(1), Fr(1, 2)} or S <= {Fr(0), Fr(1), Fr(-1)}
    raise ValueError(mode)


def nice_ok(item, ROOTS):
    f, v, rr, a0 = item
    if not path_ok(f, ROOTS, allow_above=True):
        return False
    if a0 == 'inf':
        return False
    if a0 == 'zero':
        return True
    return Fr(a0) in {Fr(0), Fr(1, 2), Fr(-1), Fr(1)}


def rho(S):
    memo = {}

    def r(word):
        if word in memo:
            return memo[word]
        n = len(word)
        if n == 1:
            res = {word: Fr(1)}
        else:
            res = defaultdict(Fr)
            for w, v in r(word[:-1]).items():
                res[w + (word[-1],)] += Fr(n - 1, n) * v
            for w, v in r(word[1:]).items():
                res[w + (word[0],)] -= Fr(n - 1, n) * v
            res = {k: v for k, v in res.items() if v}
        memo[word] = res
        return res
    out = defaultdict(Fr)
    for w, v in S.items():
        for u, x in r(w).items():
            out[u] += v * x
    return {k: v for k, v in out.items() if v}


def solve(cols, SF):
    """cols: list of ('N', term) / ('L', (f, v, rho)); SF: target symbol.  Exact solve; None if unsolvable."""
    obF, rhoF = delta22(SF), rho(SF)
    ob_keys = sorted(set(obF) | set().union(*[set(t['ob']) for k, t in cols if k == 'N']))
    labels = sorted({t['pt'] for k, t in cols if k == 'N'})
    lab_idx = {q: i for i, q in enumerate(labels)}
    rho_keys = sorted(set(rhoF) | set().union(*[set(t['rho']) for k, t in cols if k == 'N'])
                      | set().union(*[set(t[2]) for k, t in cols if k == 'L']))
    nrow = len(ob_keys) + len(labels) + len(rho_keys)
    ncol = len(cols)
    M = flint.fmpq_mat(nrow, ncol + 1)

    def setv(i, j, v):
        v = Fr(v)
        M[i, j] = flint.fmpq(v.numerator, v.denominator)
    oi = {k: i for i, k in enumerate(ob_keys)}
    ri = {k: len(ob_keys) + len(labels) + i for i, k in enumerate(rho_keys)}
    for j, (kind, t) in enumerate(cols):
        if kind == 'N':
            for k, v in t['ob'].items():
                setv(oi[k], j, v)
            setv(len(ob_keys) + lab_idx[t['pt']], j, 1)
            for k, v in t['rho'].items():
                setv(ri[k], j, v)
        else:
            for k, v in t[2].items():
                setv(ri[k], j, v)
    for k, v in obF.items():
        setv(oi[k], ncol, v)
    for k, v in rhoF.items():
        setv(ri[k], ncol, v)
    R, rank = M.rref()
    piv = []
    for i in range(rank):
        for j in range(ncol + 1):
            if R[i, j] != 0:
                piv.append(j)
                break
    if ncol in piv:
        return None
    sol = [Fr(0)] * ncol
    for i, j in enumerate(piv):
        q = R[i, ncol]
        sol[j] = Fr(int(q.p), int(q.q))
    return sol
