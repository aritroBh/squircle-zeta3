"""ir_utils.py -- integer relation finding via flint LLL, with full-precision verification."""
import mpmath as mp
from flint import fmpz_mat, fmpz

def _to_fmpz_mat(rows):
    n, m = len(rows), len(rows[0])
    A = fmpz_mat(n, m)
    for i in range(n):
        for j in range(m):
            A[i, j] = fmpz(int(rows[i][j]))
    return A

def find_relations(vals, dps_search=None, Cexp=None, maxcoeff=10**8, resid_decades=18, max_report=4):
    """Find integer relations sum m_i vals[i] = 0.
    vals: list of mpf (assumed accurate to mp.mp.dps).
    dps_search: decimal digits to use for the search truncation (default current dps - 8).
    Cexp: exponent for the scaling C = 10^Cexp (default dps_search - 5).
    A candidate row (m_1..m_n, R) is accepted if max|m_i| <= maxcoeff and
    |R| <= max(|m|) * 10^(-resid_decades) * max_i |round(C v_i)| / C ... we use:
    |residual| computed at full precision must be <= 10^(-resid_decades) * scale, scale = max|vals|.
    Returns list of (relation, maxcoeff, residual)."""
    n = len(vals)
    if dps_search is None:
        dps_search = mp.mp.dps - 8
    if Cexp is None:
        Cexp = dps_search - 5
    C = mp.mpf(10)**Cexp
    ci = [int(mp.nint(C*v)) for v in vals]
    rows = []
    for i in range(n):
        r = [0]*n + [ci[i]]
        r[i] = 1
        rows.append(r)
    A = _to_fmpz_mat(rows)
    B = A.lll(delta=0.99, eta=0.501)
    scale = max(abs(v) for v in vals)
    out = []
    for i in range(n):
        row = [int(B[i, j]) for j in range(n)]
        R = int(B[i, n])
        mc = max(abs(x) for x in row)
        if mc == 0 or mc > maxcoeff:
            continue
        res = mp.fsum(mp.mpf(x)*v for x, v in zip(row, vals))
        if abs(res) <= mp.mpf(10)**(-resid_decades) * scale:
            out.append((row, mc, res))
        if len(out) >= max_report:
            break
    return out

def prune_basis(vals, names=None, dps_search=None, maxcoeff=10**8, resid_decades=18, verbose=True):
    """Remove elements until no relation with |coeff|<=maxcoeff is found.
    Returns (kept_vals, kept_names, dependencies)."""
    if names is None:
        names = [f"v{i}" for i in range(len(vals))]
    vals = list(vals); names = list(names); deps = []
    while True:
        rels = find_relations(vals, dps_search=dps_search, maxcoeff=maxcoeff,
                              resid_decades=resid_decades, max_report=1)
        if not rels:
            break
        row, mc, res = rels[0]
        nz = [i for i, x in enumerate(row) if x]
        j = nz[-1]  # eliminate last involved element
        deps.append((names[j], list(zip([names[i] for i in nz], [row[i] for i in nz])), mc, res))
        if verbose:
            print(f"    prune: remove {names[j]}  (relation coeffs {[(names[i], row[i]) for i in nz]}, maxc {mc}, resid {mp.nstr(res,2)})")
        del vals[j]; del names[j]
    return vals, names, deps
