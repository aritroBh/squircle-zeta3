"""Group g6, common code.  Exact formal expansion of the flat-angle coefficients mu_{k,n}.

Written from the formulas of the text (sec_g6_flat.tex), not from any file of the notes.

  mu_{k,n} = -1/(n(2n-1)) * sum_{j,i,m} e_{j,i,m}(n) * (n)_j/(2n)_j * D_{i,m}(n+j,n),

where [u^k]E_n = sum e_{j,i,m}(n) y^j a^i b^m,  E_n = exp(uA) (1+y(exp(u l)-1))^(1-2n),
A = -(n-1)a + n b, l = b - a, and D_{i,m}(p,q) = i! m! [x^i x'^m] exp(Delta),
Delta = Delta^H + Delta^Z (zeta values are formal variables Z_w).
Returned objects: dict {exponent tuple (e_2,...,e_K): Fraction}; the monomial is prod Z_w^{e_w}.
"""
from fractions import Fraction
from math import factorial
from sympy import QQ, ring


def harm(r, N):
    """H^{(r)}_N = sum_{i=1}^N i^{-r} as a Fraction."""
    return sum((Fraction(1, i ** r) for i in range(1, N + 1)), Fraction(0))


def E_coeffs(n, K):
    """List over k=0..K of dict {(j,i,m): Fraction} = coefficients of y^j a^i b^m in [u^k] E_n."""
    R, u, y, a, b = ring("u,y,a,b", QQ)

    def trunc(p):
        return R.from_dict({m: c for m, c in p.items() if m[0] <= K})

    ell = b - a
    A = -(n - 1) * a + n * b
    # X = y (exp(u ell) - 1)
    X = R.zero
    term = R.one
    for r in range(1, K + 1):
        term = trunc(term * u * ell) * QQ(1, r)
        X = X + y * term
    L = R.zero
    Xp = R.one
    for m in range(1, K + 1):
        Xp = trunc(Xp * X)
        L = L + Xp * QQ((-1) ** (m - 1), m)
    logE = u * A + (1 - 2 * n) * L
    E = R.one
    P = R.one
    for t in range(1, K + 1):
        P = trunc(P * logE) * QQ(1, t)
        E = E + P
    out = [dict() for _ in range(K + 1)]
    for (du, dy, da, db), c in E.items():
        out[du][(dy, da, db)] = Fraction(int(c.numerator), int(c.denominator))
    return out


def exp_delta(p, q, K):
    """exp(Delta(p,q;x,x')) truncated at total degree K, as dict {(i,m,(e_2..e_K)): Fraction}."""
    names = ["x", "xp"] + ["Z%d" % w for w in range(2, K + 1)]
    R = ring(",".join(names), QQ)
    Rr, gens = R[0], R[1:]
    x, xp = gens[0], gens[1]
    Z = {w: gens[w] for w in range(2, K + 1)}  # Z[w] = gens[w] since gens[2] is Z2

    def deg(m):
        return m[0] + m[1]

    def trunc(P):
        return Rr.from_dict({m: c for m, c in P.items() if deg(m) <= K})

    Delta = Rr.zero
    s = x + xp
    for r in range(1, K + 1):
        sgn = QQ((-1) ** (r - 1), r)
        hp = QQ(*_fr(harm(r, p - 1)))
        hq = QQ(*_fr(harm(r, q - 1)))
        hpq = QQ(*_fr(harm(r, p + q - 1)))
        Delta = Delta + sgn * (hp * x ** r + hq * xp ** r - hpq * s ** r)
    for w in range(2, K + 1):
        Delta = Delta + QQ((-1) ** (w + 1), w) * Z[w] * (s ** w - x ** w - xp ** w)
    E = Rr.one
    P = Rr.one
    for t in range(1, K + 1):
        P = trunc(P * Delta) * QQ(1, t)
        E = E + P
    out = {}
    for m, c in E.items():
        out[(m[0], m[1], tuple(m[2:]))] = Fraction(int(c.numerator), int(c.denominator))
    return out


def _fr(f):
    return (f.numerator, f.denominator)


def poch_ratio(n, j):
    """(n)_j / (2n)_j"""
    r = Fraction(1)
    for t in range(j):
        r *= Fraction(n + t, 2 * n + t)
    return r


def mu_formal(n, K):
    """List over k=0..K of dict {zeta exponent tuple (e_2..e_K): Fraction}, the formal mu_{k,n}."""
    Ec = E_coeffs(n, K)
    EX = {j: exp_delta(n + j, n, K) for j in range(K + 1)}
    pref = Fraction(-1, n * (2 * n - 1))
    res = []
    for k in range(K + 1):
        acc = {}
        for (j, i, m), c in Ec[k].items():
            r = poch_ratio(n, j)
            fac = Fraction(factorial(i) * factorial(m))
            for (ii, mm, ze), v in EX[j].items():
                if ii == i and mm == m:
                    acc[ze] = acc.get(ze, 0) + pref * c * r * fac * v
        res.append({z: v for z, v in acc.items() if v != 0})
    return res


def zeros(K):
    return tuple([0] * (K - 1))


def mono(K, **kw):
    """exponent tuple from keywords like z2=1,z3=2"""
    e = [0] * (K - 1)
    for key, v in kw.items():
        e[int(key[1:]) - 2] = v
    return tuple(e)


if __name__ == "__main__":
    K = 3
    for n in (1, 2, 3):
        print(n, mu_formal(n, K))
