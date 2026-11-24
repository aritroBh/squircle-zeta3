"""LLL integer-relation helper (python-flint): returns the short vectors whose scaled dot product with the numbers is ~0."""
import mpmath as mp
from flint import fmpz_mat
def lll_relations(vals, digits, maxshow=12):
    """vals: list of mpf (current mp.mp.dps >= digits+10). Returns list of (relation(list of int), residual mpf), best first."""
    n = len(vals)
    scale = mp.mpf(10) ** digits
    rows = []
    for i, x in enumerate(vals):
        row = [0] * n + [int(mp.nint(x * scale))]
        row[i] = 1
        rows.append(row)
    M = fmpz_mat(rows).lll()
    out = []
    for r in range(M.nrows()):
        rel = [int(M[r, c]) for c in range(n)]
        last = int(M[r, n])
        res = sum(a * x for a, x in zip(rel, vals))
        out.append((rel, res, last))
    return out

def find_rel(target, basis, D, maxshow=3, coeffcap=None):
    """LLL search for relations among [target]+[basis values] (basis: list of (name, mpf)); mp.mp.dps must be >= D+15.
    Returns list of dicts for every reduced-basis row whose residual is < 10^-(D-25); sorted by max|coeff|."""
    vals = [target] + [v for _, v in basis]
    names = ["TARGET"] + [n for n, _ in basis]
    rows = lll_relations(vals, D)
    out = []
    for rel, res, last in rows:
        if abs(res) < mp.mpf(10) ** (-(D - 25)):
            out.append({"rel": rel, "res": res, "maxcoeff": max(abs(a) for a in rel),
                        "terms": [(names[i], a) for i, a in enumerate(rel) if a != 0]})
    out.sort(key=lambda d: d["maxcoeff"])
    return out


def find_rel2(target, basis, D):
    """like find_rel but separates relations with a nonzero target coefficient from internal ones (target coefficient 0)."""
    R = find_rel(target, basis, D)
    return [r for r in R if r["rel"][0] != 0], [r for r in R if r["rel"][0] == 0]
