code = '''"""Exact linear algebra for the HEXC2 certificate: select independent proved rows mod p,
solve y . S = G exactly over Q, verify exactly."""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import nmod_mat, fmpz_mat, fmpq_mat, ctx
ctx.threads = 1

PRIMES = [1000000009, 1000000097, 1000000241]

def cmod(c, p):
    return (c.numerator % p) * pow(c.denominator % p, -1, p) % p if hasattr(c, 'numerator') else c % p

def rows_to_nmod(rows, ncols, ids, p):
    mat = nmod_mat(len(rows), ncols, p)
    for i, row in enumerate(rows):
        for w, c in row:
            j = ids.get(w)
            if j is not None: mat[i, j] = cmod(c, p)
    return mat

def pivots_of_rref(rr, rank):
    """pivot column indices of an RREF matrix (first `rank` rows nonzero)."""
    pp = []
    for i in range(rank):
        row = rr[i]
        for j in range(rr.ncols()):
            if row[j] != 0:
                pp.append(j); break
    return pp

def select_independent_rows(rows, cols, p=PRIMES[0]):
    """indices of original rows forming a basis of the row space mod p: pivot COLUMNS of rref(M^T)."""
    ids = {w: i for i, w in enumerate(cols)}
    MT = nmod_mat(len(cols), len(rows), p)
    for i, row in enumerate(rows):
        for w, c in row:
            j = ids.get(w)
            if j is not None: MT[j, i] = cmod(c, p)
    rr, r = MT.rref(inplace=True)
    prow = pivots_of_rref(rr, r)
    # sanity: the selected rows are independent
    S = rows_to_nmod([rows[i] for i in prow], len(cols), ids, p)
    assert S.rank() == r, "selected rows not independent"
    return prow, r

def select_pivot_cols(sel_rows, cols, p=PRIMES[0]):
    ids = {w: i for i, w in enumerate(cols)}
    S = rows_to_nmod(sel_rows, len(cols), ids, p)
    rr, r = S.rref(inplace=True)
    return pivots_of_rref(rr, r)

def exact_solve(sel_rows, pcols, cols, target):
    """Solve y . S = G exactly. S = sel_rows (r x m), G = target dict; uses pivot columns."""
    r = len(sel_rows); m = len(cols)
    ids = {w: i for i, w in enumerate(cols)}
    # Build B = S[:, pcols] as fmpz (rows scaled to integers) and Gv = target[pcols]
    from math import lcm
    B = []   # r x r exact Fractions
    dB = []
    for row in sel_rows:
        d = {w: c for w, c in row}
        B.append([d.get(cols[j], Fr(0)) for j in pcols])
    den = 1
    for brow in B:
        for c in brow: den = lcm(den, c.denominator)
    Bz = fmpz_mat([[int(c * den) for c in brow] for brow in B])
    Gv = [target.get(cols[j], Fr(0)) for j in pcols]
    denG = 1
    for c in Gv: denG = lcm(denG, c.denominator)
    bz = fmpz_mat([[int(c * denG)] for c in Gv])   # r x 1
    # solve Bz^T . x = (den/denG?) -- careful: y . B = Gv  =>  B^T y^T = Gv^T
    # with integer matrices: (den B)^{-T}: solve Bz^T X = bz, then y = (denG/den) X
    sol, sden = fmpz_mat.solve(Bz.transpose(), bz)
    # y = (denG/den) * sol / sden  (as rationals); verify
    y = [Fr(int(sol[i, 0]) * denG, int(sden) * den) for i in range(r)]
    return y

def verify_combination(rows, ys, cols, target):
    """exact check: sum_i y_i rows[i] == target on all columns."""
    acc = {}
    for row, y in zip(rows, ys):
        if y == 0: continue
        for w, c in row:
            acc[w] = acc.get(w, Fr(0)) + c * y
    for w in cols:
        if acc.get(w, Fr(0)) != target.get(w, Fr(0)):
            return False, w
    extra = set(target) - set(cols)
    return (len(extra) == 0), None if not extra else ('target words outside cols', extra)
'''
open('exact_linalg.py','w').write(code)
print('written exact_linalg.py')
