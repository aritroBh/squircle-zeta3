"""g2_fit.py -- shared machinery for the B4 HPL fits (Task 2).

HPLs H_w(c), letters {0,1,-1}, H_{a,w}(c) = int_0^c f_a H_w, f_0 = 1/t, f_1 = 1/(1-t), f_{-1} = 1/(1+t).
hpl_series(w, order): exact Taylor coefficients [c^m] H_w(c).
Even-sector reduction: the target M4(c) = sum mu_{4,n} c^{2n} is even. For each orbit {w,-w} and
exponent e, the even part of c^e H_w is proportional to
  F_{w,e}(c) = c^e ( H_w(c) + (-1)^{e+nz(w)} H_{-w}(c) ),   nz(w) = #nonzero letters of w.
A solution exists in the full basis iff one exists in the even basis (projection argument).
solve_exact: exact fmpq rref on even Taylor rows; accept only if every equation is satisfied.
"""
import itertools, json
from fractions import Fraction as Fr
import flint


def words(maxlen):
    W = []
    for L in range(1, maxlen + 1):
        for w in itertools.product([0, 1, -1], repeat=L):
            if w[-1] != 0:
                W.append(w)
    return W


def nz(w):
    return sum(1 for a in w if a != 0)


def neg(w):
    return tuple(-a for a in w)


def orbit_rep(w):
    return min(w, neg(w))


from functools import lru_cache


@lru_cache(maxsize=None)
def hpl_series(w, order):
    """[c^m] H_w(c) for m = 0..order (exact)."""
    w = tuple(w)
    N = order + 3
    if len(w) == 1:
        a = w[0]
        s = [Fr(0)]*(N + 1)
        for m in range(1, N + 1):
            s[m] = Fr(1, m) if a == 1 else Fr((-1)**(m + 1), m)
        return s
    inner = hpl_series(tuple(w[1:]), order)
    a = w[0]
    out = [Fr(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m]/m
    else:
        acc = Fr(0)
        for m in range(0, N):
            acc = acc*a + inner[m]
            out[m + 1] = acc/(m + 1)
    return out


def even_series(w, e, order):
    """coefficients of c^e (H_w + (-1)^{e+nz(w)} H_{-w}) as list [0..order]."""
    s1 = hpl_series(w, order)
    w2 = neg(w)
    s2 = hpl_series(w2, order) if w2 != w else s1
    sgn = 1 if (e + nz(w)) % 2 == 0 else -1
    out = [Fr(0)]*(order + 1)
    for j in range(max(0, e), order + 1):
        m = j - e                       # m in [0, order - e]; s1,s2 have length order + 4
        out[j] = s1[m] + sgn*s2[m]
    return out


def build_even_basis(maxlen, exps, order, wfilter=None):
    """list of (name, series) for the even-sector basis, plus the constant '1'."""
    seen = set()
    basis = [('1', [Fr(1)] + [Fr(0)]*order)]
    for w in words(maxlen):
        if wfilter and not wfilter(w):
            continue
        r = orbit_rep(w)
        if r in seen:
            continue
        seen.add(r)
        for e in exps:
            ser = even_series(w, e, order)
            basis.append((f'c^{e} He{list(w)}', ser))
    return basis


def target_series(part, order, tablefile='flat_exact_N520_K4_formula.json'):
    data = json.load(open(tablefile))['4']
    T = [Fr(0)]*(order + 1)
    for nn, P in data.items():
        if 2*int(nn) <= order:
            T[2*int(nn)] = Fr(P.get(part, '0'))
    return T


def solve_exact(basis, T, order, verbose=True):
    """exact solve on the even rows 0,2,...,order (and odd rows must be 0: they are, by construction)."""
    rows = list(range(0, order + 1, 2))
    nb = len(basis)
    if verbose:
        print(f'   {nb} basis functions, {len(rows)} even equations')
    Mx = flint.fmpq_mat(len(rows), nb)
    for j, (_, ser) in enumerate(basis):
        for i, r in enumerate(rows):
            v = ser[r]
            Mx[i, j] = flint.fmpq(v.numerator, v.denominator)
    rhs = flint.fmpq_mat(len(rows), 1)
    for i, r in enumerate(rows):
        rhs[i, 0] = flint.fmpq(T[r].numerator, T[r].denominator)
    aug = flint.fmpq_mat(len(rows), nb + 1)
    for i in range(len(rows)):
        for j in range(nb):
            aug[i, j] = Mx[i, j]
        aug[i, nb] = rhs[i, 0]
    R, rank = aug.rref()
    rankA = Mx.rref()[1]
    if verbose:
        print('   rank(A) =', rankA, ' rank([A|b]) =', rank)
    if rank != rankA:
        return None, rankA, rank
    sol = [flint.fmpq(0)]*nb
    for i in range(rank):
        row = [R[i, j] for j in range(nb + 1)]
        piv = next(j for j in range(nb + 1) if row[j] != 0)
        if piv == nb:
            return None, rankA, rank
        sol[piv] = row[nb]
    terms = [(basis[j][0], sol[j]) for j in range(nb) if sol[j] != 0]
    return terms, rankA, rank
