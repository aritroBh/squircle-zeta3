sunit = r'''"""cor6_sunit.py -- S-unit bases via norm certification + log lattice + exact verification."""
from fractions import Fraction as Fr
import mpmath as mp
from itertools import product as iproduct

def smooth_part(q, S0):
    """q = Fr; if q = prod p^e over S0 return exponent dict; else None."""
    e = {}
    for p in S0:
        a = 0
        n, d = q.numerator, q.denominator
        while n % p == 0: n //= p; a += 1
        while d % p == 0: d //= p; a -= 1
        if a: e[p] = a
    n, d = q.numerator, q.denominator
    for p in S0:
        while n % p == 0: n //= p
        while d % p == 0: d //= p
    if n == 1 and d == 1:
        return e
    return None

def is_sunit(F, a, S0):
    """exact certificate: a is an S0-unit iff N(a) is a S0-number (times sign)."""
    N = F.norm(a)
    if N == 0: return None
    return smooth_part(abs(N), S0) is not None

def find_sunits(F, S0, bound=3, denom=2, want_units_only=False, limit=4000):
    """search elements a0+a1 t+... with ai in (1/denom)Z, |ai|<=bound, S-units."""
    out = []
    d = F.d
    rng = [Fr(k, denom) for k in range(-bound * denom, bound * denom + 1)]
    cnt = 0
    for coeffs in iproduct(rng, repeat=d):
        if not any(coeffs): continue
        a = tuple(coeffs)
        N = F.norm(a)
        if N == 0: continue
        if want_units_only and abs(N) != 1: continue
        if smooth_part(abs(N), S0) is not None:
            out.append(a)
            cnt += 1
            if cnt >= limit: break
    return out

def log_matrix(F, elems, dps=120):
    rts = F.roots(dps)
    M = []
    for a in elems:
        row = []
        for r in rts:
            v = F.eval(a, r)
            row.append(mp.log(abs(v)))
        M.append(row)
    return M

def num_rank(M, tol=None):
    mp.mp.dps = 100
    A = mp.matrix(M)
    if tol is None:
        tol = mp.mpf(10) ** (-60)
    # Gaussian elimination
    rows = A.rows; cols = A.cols
    R = 0
    for c in range(cols):
        p = max(range(R, rows), key=lambda r: abs(A[r, c]))
        if abs(A[p, c]) < tol: continue
        A.swap_row(R, p)
        for r in range(R + 1, rows):
            if A[r, c]:
                A[r, :] -= (A[r, c] / A[R, c]) * A[R, :]
        R += 1
        if R == rows: break
    return R

def exponent_vector(F, x, basis, dps=150, max_denom=10**6):
    """solve Lambda(x) = sum e_i Lambda(g_i); rational-reconstruct; EXACT verify.
    Returns exponent tuple (Fr), or None if failed."""
    rts = F.roots(dps)
    bx = [mp.log(abs(F.eval(x, r))) for r in rts]
    B = [[mp.log(abs(F.eval(g, r))) for g in basis] for r in rts]
    # least squares via normal equations (well-conditioned for small r)
    r = len(basis)
    G = [[sum(B[k][i] * B[k][j] for k in range(len(rts))) for j in range(r)] for i in range(r)]
    b = [sum(B[k][i] * bx[k] for k in range(len(rts))) for i in range(r)]
    # solve
    A = [row[:] + [b[i]] for i, row in enumerate(G)]
    for col in range(r):
        p = max(range(col, r), key=lambda rr: abs(A[rr][col]))
        A[col], A[p] = A[p], A[col]
        for rr in range(col + 1, r):
            A[rr] = [x - y * A[rr][col] / A[col][col] for x, y in zip(A[rr], A[col])]
    e = [mp.mpf(0)] * r
    for i in range(r - 1, -1, -1):
        e[i] = (A[i][r] - sum(A[i][j] * e[j] for j in range(i + 1, r))) / A[i][i]
    eq = []
    for v in e:
        f = Fr(v).limit_denominator(max_denom)
        if abs(mp.mpf(f.numerator) / f.denominator - v) > mp.mpf(10) ** (-30):
            return None
        eq.append(f)
    # exact verify: x == +- prod g_i^e_i  (clear denominators of exponents)
    from math import gcd
    L = 1
    for f in eq: L = L * f.denominator // gcd(L, f.denominator)
    # check x^L == +- prod g_i^(e_i*L)
    lhs = F.pow(x, L)
    rhs = F.one()
    for f, g in zip(eq, basis):
        rhs = F.mul(rhs, F.pow(g, int(f * L)))
    if F.eq(lhs, rhs) or F.eq(lhs, F.neg(rhs)):
        return tuple(eq)
    return None

def build_basis(F, S0, dimV, extra=(), search_bound=3, search_denom=2, dps=150, verbose=True):
    """find dimV certified S-units, numerically independent, spanning V (log rank = dimV),
    including exact exponent vectors for the 'extra' required elements."""
    pool = list(extra)
    # always include 2 and small rational primes
    for p in S0:
        pool.append(F.from_int(p))
    found = find_sunits(F, S0, bound=search_bound, denom=search_denom, limit=6000)
    pool += found
    # greedy independent set by log rank
    basis = []
    cur = []
    for a in pool:
        M = log_matrix(F, basis + [a], dps=80)
        if num_rank(M) > len(basis):
            basis.append(a)
            if len(basis) == dimV: break
    if verbose:
        print(f"  basis size {len(basis)} (dim V target {dimV})")
    if len(basis) < dimV:
        return None, None
    # check no small relations: exponent_vector of each basis element should be trivial unit vector
    return basis, pool
'''
with open("cor6_sunit.py", "w") as f:
    f.write(sunit)
print("ok")
