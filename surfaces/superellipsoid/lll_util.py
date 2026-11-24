# lll_util.py -- integer relation search by LLL (python-flint), as in main.tex.
import mpmath as mp
import flint

def lll_relation(vec, digits, check_digits=None):
    """vec: list of mpf (target first). Returns (relation list, residual) of the shortest LLL vector."""
    n = len(vec)
    S = mp.mpf(10)**digits
    rows = []
    for i, x in enumerate(vec):
        row = [0]*n
        row[i] = 1
        row.append(int(mp.nint(x*S)))
        rows.append(row)
    M = flint.fmpz_mat(rows)
    R = M.lll()
    best = None
    for i in range(R.nrows()):
        rel = [int(R[i, j]) for j in range(n)]
        if rel[0] == 0: continue
        res = mp.fsum(c*x for c, x in zip(rel, vec))
        h = max(abs(c) for c in rel)
        if best is None or h < best[2]:
            best = (rel, res, h)
    return best

def describe(rel, names):
    c0 = rel[0]
    return "(1/%d)*[%s]" % (-c0, " + ".join("%d*%s" % (c, n) for c, n in zip(rel[1:], names) if c))
