"""m5fit.py -- FIT machinery (not part of any proof).

Finds rational coefficients q_{e,w} with
    sum_n T_n c^{2n} = q_0 + q_{-2} c^{-2} + sum_{e,w} q_{e,w} c^e (H_w(c) + s H_{-w}(c)),   s = (-1)^{e+nz(w)},
for a target sequence T_n (n >= 1) given by a Z-sum dict {word: Q(n)} (from mu5_zsum.pkl).
Each term c^e(H_w(c)+sH_{-w}(c)) = c^e H_w(c) + (-c)^e H_w(-c) is even in c, so only even coefficients are fitted:
    [c^{2n}] = 2 a_w[2n-e]  (a_w = Taylor coefficients of H_w),  n = 1..Nn,   and [c^0] = [c^{-2}] = 0.
Words: last letter 1 (representative of the pair {w,-w}), any letters 0,+-1 before, length 1..W.
All arithmetic mod a prime p (flint nmod_mat); solutions are lifted by CRT + rational reconstruction over
several primes and must agree.  The exact proof is a separate script (s03_prove.py).
"""
import pickle, itertools, sys
from fractions import Fraction as Fr
import sympy as sp
import flint

PRIMES = [2305843009213693951, 4611686018427387847, 9223372036854775783, 4611686018427387817]


def words_upto(W):
    out = []
    for L in range(1, W + 1):
        for pre in itertools.product((0, 1, -1), repeat=L - 1):
            out.append(tuple(pre) + (1,))
    return out


def load_mu(fname="mu5_zsum.pkl"):
    n = sp.Symbol('n')
    raw = pickle.load(open(fname, "rb"))
    out = {}
    for zm, E in raw.items():
        out[zm] = {w: (sp.Poly(sp.sympify(a), n, domain='QQ'), sp.Poly(sp.sympify(b), n, domain='QQ')) for w, (a, b) in E.items()}
    return out


def poly_mod(P, x, p):
    v = 0
    for c in P.all_coeffs():
        c = sp.Rational(c)
        v = (v * x + int(c.p) * pow(int(c.q), -1, p)) % p
    return v


def Zvals(word, Mmax, p, inv, cache):
    """list Z_word(M) mod p, M = 0..Mmax."""
    if word in cache:
        return cache[word]
    if not word:
        r = [1] * (Mmax + 1)
    else:
        a, rest = word[0], word[1:]
        R = Zvals(rest, Mmax, p, inv, cache)
        r = [0] * (Mmax + 1)
        ea = abs(a)
        for M in range(1, Mmax + 1):
            t = pow(inv[M], ea, p) * R[M - 1] % p
            if a < 0 and M % 2:
                t = -t
            r[M] = (r[M - 1] + t) % p
    cache[word] = r
    return r


def target_mod(E, Nn, p, inv):
    """T_n mod p for n = 1..Nn from the Z-dict E (values of polys)."""
    cache = {}
    T = [0] * (Nn + 1)
    for w, (num, den) in E.items():
        Z = Zvals(w, 2 * Nn, p, inv, cache)
        for n in range(1, Nn + 1):
            T[n] = (T[n] + poly_mod(num, n, p) * pow(poly_mod(den, n, p), -1, p) % p * Z[2 * n]) % p
    return T


def hpl_coeffs_mod(words, Nmax, p, inv):
    """a_w[m] mod p, m=0..Nmax, for every word in `words` (and its suffixes)."""
    memo = {}
    def get(w):
        if w in memo:
            return memo[w]
        if len(w) == 1:
            a = w[0]
            r = [0] + [inv[m] if a == 1 else (inv[m] if m % 2 else p - inv[m]) for m in range(1, Nmax + 1)]
        else:
            inner = get(w[1:]); a = w[0]
            r = [0] * (Nmax + 1)
            if a == 0:
                for m in range(1, Nmax + 1):
                    r[m] = inner[m] * inv[m] % p
            else:
                acc = 0
                for m in range(1, Nmax + 1):
                    acc = (acc * a + inner[m - 1]) % p
                    r[m] = acc * inv[m] % p
        memo[w] = r
        return r
    return {w: get(w) for w in words}


def build_system(E, W, erange, Nn, p):
    inv = [0] + [pow(m, -1, p) for m in range(1, 2 * Nn + 10)]
    words = words_upto(W)
    cols = [('const', 0), ('const', -2)] + [('He', e, w) for e in erange for w in words]
    A = hpl_coeffs_mod(words, 2 * Nn + 6, p, inv)
    T = target_mod(E, Nn, p, inv)
    rows = []
    rhs = []
    for N in [-2, 0] + [2 * n for n in range(1, Nn + 1)]:
        row = []
        for col in cols:
            if col[0] == 'const':
                row.append(1 if N == col[1] else 0)
            else:
                _, e, w = col
                m = N - e
                row.append(2 * A[w][m] % p if m >= 1 else 0)
        rows.append(row)
        rhs.append(0 if N <= 0 else T[N // 2])
    return cols, rows, rhs


def solve_mod(rows, rhs, p):
    nr, nc = len(rows), len(rows[0])
    M = flint.nmod_mat(nr, nc + 1, [x for r, b in zip(rows, rhs) for x in (r + [b])], p)
    R, rank = M.rref()
    Rl = R.tolist()
    piv = []
    for i in range(rank):
        row = Rl[i]
        j = next(j for j in range(nc + 1) if int(row[j]) != 0)
        piv.append(j)
    if nc in piv:
        return None, rank, None
    sol = [0] * nc
    for i, j in enumerate(piv):
        sol[j] = int(Rl[i][nc])
    return sol, rank, piv


def ratrec(a, m):
    """rational reconstruction of a mod m."""
    import math
    bound = math.isqrt(m // 2)
    r0, r1, s0, s1 = m, a % m, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def fit(E, W, erange, Nn, nprimes=3, verbose=True):
    sols = []; mods = []
    info = None
    for p in PRIMES[:nprimes]:
        cols, rows, rhs = build_system(E, W, erange, Nn, p)
        sol, rank, piv = solve_mod(rows, rhs, p)
        if verbose:
            print(f"   p={p}: {len(rows)} equations, {len(cols)} unknowns, rank {rank}, consistent {sol is not None}", flush=True)
        if sol is None:
            return None, cols, rank
        info = (rank, len(cols), piv)
        sols.append(sol); mods.append(p)
    # CRT
    from sympy.ntheory.modular import crt
    out = []
    m = 1
    for x in mods: m *= x
    for j in range(len(cols)):
        v, _ = crt(mods, [s[j] for s in sols])
        out.append(ratrec(int(v), m))
    # consistency: reconstruct from first nprimes-1 primes too
    m2 = 1
    for x in mods[:-1]: m2 *= x
    out2 = []
    for j in range(len(cols)):
        v, _ = crt(mods[:-1], [s[j] for s in sols[:-1]])
        out2.append(ratrec(int(v), m2))
    stable = out == out2
    if verbose:
        print(f"   rational reconstruction stable over {len(mods)-1} vs {len(mods)} primes: {stable}; nonzero {sum(1 for x in out if x)}")
    return {cols[j]: out[j] for j in range(len(cols)) if out[j]}, cols, info


def to_json_terms(sol):
    out = []
    for col, q in sol.items():
        if col[0] == 'const':
            out.append(['1' if col[1] == 0 else f'c^{col[1]}', str(q)])
        else:
            _, e, w = col
            out.append([f"c^{e} He[{','.join(str(a) for a in w)}]", str(q)])
    return out
