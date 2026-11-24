
cert_module = '''"""thx1_cert.py -- certificates of MZV identities from proved relation instances.

Each relation instance R_j (DS: finite double shuffle from absolute convergence;
DUAL: change of variables t->(1-t)/(1+t)) is a PROVED identity among word-values.
A certificate for target V: exact rationals h_j with V = sum_j h_j R_j (exact sparse
Fraction verification) => value(V) = 0.

Structure via flint mod p: rref(R) -> pivot columns; rref(R^T) -> independent original rows.
Solve h mod primes (square pivot system, batch), CRT + rational reconstruction, exact verify.
"""
from fractions import Fraction as Fr
import flint, pickle, os

def build_nmod(rels, cols, p):
    colpos = {c: i for i, c in enumerate(cols)}
    A = flint.nmod_mat(len(rels), len(cols), p)
    for i, r in enumerate(rels):
        for k, v in r.items():
            A[i, colpos[k]] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    return A, colpos

def rref_pivots(A):
    R = A.rref()
    piv = []
    for i in range(R.nrows()):
        row_lead = None
        for j in range(R.ncols()):
            if R[i, j] != 0:
                row_lead = j
                break
        if row_lead is not None:
            piv.append(row_lead)
    return piv, R

def compute_structure(tag, rels, cols, p, recompute=False):
    fn = f'struct_{tag}.pkl'
    if os.path.exists(fn) and not recompute:
        return pickle.load(open(fn, 'rb'))
    A, colpos = build_nmod(rels, cols, p)
    pivcols, _ = rref_pivots(A)
    At = A.transpose()
    pivrows, _ = rref_pivots(At)
    assert len(pivcols) == len(pivrows)
    st = {'rank': len(pivcols), 'pivcols': pivcols, 'indep_rows': pivrows, 'colpos': colpos, 'prime': p}
    pickle.dump(st, open(fn, 'wb'))
    return st

def solve_h_mod_batch(targets, rels, st, cols, p):
    """targets: list of exact Fraction combos. Returns list of h (len r each) mod p."""
    B = st['indep_rows']; pivcols = st['pivcols']; colpos = st['colpos']
    r = len(B)
    posmap = {c: j for j, c in enumerate(pivcols)}
    S = flint.nmod_mat(r, r, p)
    for i, bi in enumerate(B):
        for k, v in rels[bi].items():
            j = posmap.get(colpos[k])
            if j is not None:
                S[j, i] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    b = flint.nmod_mat(r, len(targets), p)
    for tj, V in enumerate(targets):
        for k, v in V.items():
            j = posmap.get(colpos[k])
            if j is not None:
                b[j, tj] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    x = S.solve(b)
    return [[int(x[i, tj]) % p for i in range(r)] for tj in range(len(targets))]

def check_h_mod(h, V, rels, st, cols, p):
    """verify sum_i h_i R_{B[i]} = V over ALL columns mod p (sparse)."""
    B = st['indep_rows']; colpos = st['colpos']
    acc = {}
    for i, hi in enumerate(h):
        if hi == 0: continue
        for k, v in rels[B[i]].items():
            vv = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
            acc[k] = (acc.get(k, 0) + hi * vv) % p
    for k, v in V.items():
        vv = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
        acc[k] = (acc.get(k, 0) - vv) % p
    return all(x % p == 0 for x in acc.values())

def crt(vals):
    """vals: list of (residue, prime). Returns (x mod P, P)."""
    P = 1
    for _, p in vals: P *= p
    x = 0
    for a, p in vals:
        q = P // p
        x = (x + a * q * pow(q, -1, p)) % P
    return x, P

def rat_recon(x, P):
    """rational reconstruction: find a/b with x = a/b mod P, |a|,b <= sqrt(P/2)."""
    if x == 0: return Fr(0)
    r0, r1 = P, x % P
    t0, t1 = 0, 1
    bound = 1
    while r1 * r1 * 2 > P:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    if t1 < 0: r1, t1 = -r1, -t1
    f = Fr(r1, t1)
    if f.denominator * f.denominator * 2 > P or abs(f.numerator) * abs(f.numerator) * 2 > P:
        return None
    return f

def exact_verify(h_rat, V, rels, B):
    """THE certificate check: sum_i h_i R_{B[i]} - V == 0 exactly (Fractions)."""
    acc = {}
    for i, hi in enumerate(h_rat):
        if hi == 0: continue
        for k, v in rels[B[i]].items():
            acc[k] = acc.get(k, Fr(0)) + hi * v
    for k, v in V.items():
        acc[k] = acc.get(k, Fr(0)) - v
    acc = {k: v for k, v in acc.items() if v}
    return acc == {}, acc

def certify(tag, V, rels, st, cols, primes, max_primes=None):
    """full pipeline for one target V. Returns (h_rat, ok) or (None, False)."""
    B = st['indep_rows']
    hmods = []
    for p in primes:
        hs = solve_h_mod_batch([V], rels, st, cols, p)
        h = hs[0]
        assert check_h_mod(h, V, rels, st, cols, p), f"mod-p check failed {tag} p={p}"
        hmods.append(h)
    # try reconstruction with growing prime count
    for k in range(2, len(primes)+1):
        Ps = primes[:k]
        h_rat = []
        ok = True
        for i in range(len(B)):
            x, P = crt([(hmods[j][i], Ps[j]) for j in range(k)])
            f = rat_recon(x, P)
            if f is None:
                ok = False; break
            h_rat.append(f)
        if not ok: continue
        good, acc = exact_verify(h_rat, V, rels, B)
        if good:
            return h_rat, True
    return None, False
'''
open('thx1_cert.py','w').write(cert_module)
print("written")
