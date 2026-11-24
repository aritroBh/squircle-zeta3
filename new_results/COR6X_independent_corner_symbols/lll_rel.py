
# lll_rel.py -- integer relation finding via flint fmpz_mat.lll (fast, exact)
import mpmath as mp
import flint

def find_relations(xs, prec_digits, max_report=6, scale_pow=None):
    n = len(xs)
    if scale_pow is None:
        scale_pow = prec_digits
    scale = mp.mpf(10) ** scale_pow
    B = []
    for i in range(n):
        row = [0] * (n + 1)
        row[i] = 1
        row[n] = int(mp.nint(xs[i] * scale))
        B.append(row)
    M = flint.fmpz_mat(B)
    R = M.lll()
    out = []
    rows = R.tolist()
    for row in rows[:max_report]:
        a = [int(c) for c in row[:n]]
        resid = mp.fsum(mp.mpf(int(c)) * x for c, x in zip(a, xs))
        out.append((a, resid))
    return out
