"""cor6_nf.py -- exact number-field arithmetic + S-unit machinery for COR6.

Elements of F = Q(theta) are tuples of Fractions (coeffs of poly in theta, deg < d).
theta: a root of minpoly m (tuple of Fractions, ascending, m[-1]=1 monic).
All arithmetic exact. Embeddings via mpmath roots (for log-lattice).
"""
from fractions import Fraction as Fr
import mpmath as mp

class NF:
    def __init__(self, minpoly, theta_approx=None, name="F"):
        """minpoly: list of Fr ascending, monic."""
        self.m = [Fr(x) for x in minpoly]
        assert self.m[-1] == 1
        self.d = len(self.m) - 1
        self.name = name
        self.theta_approx = theta_approx
    def zero(self): return (Fr(0),) * self.d
    def one(self): return (Fr(1),) + (Fr(0),) * (self.d - 1)
    def theta(self): return (Fr(0), Fr(1)) + (Fr(0),) * (self.d - 2) if self.d >= 2 else None
    def from_int(self, k):
        return (Fr(k),) + (Fr(0),) * (self.d - 1)
    def from_frac(self, q):
        return (Fr(q),) + (Fr(0),) * (self.d - 1)
    def add(self, a, b): return tuple(x + y for x, y in zip(a, b))
    def sub(self, a, b): return tuple(x - y for x, y in zip(a, b))
    def neg(self, a): return tuple(-x for x in a)
    def scale(self, q, a): return tuple(Fr(q) * x for x in a)
    def mul(self, a, b):
        d = self.d
        tmp = [Fr(0)] * (2 * d - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj: tmp[i + j] += ai * bj
        # reduce mod m
        m = self.m
        for k in range(2 * d - 2, d - 1, -1):
            c = tmp[k]
            if c:
                for i in range(d):
                    tmp[k - d + i] -= c * m[i]
        return tuple(tmp[:d])
    def pow(self, a, n):
        if n < 0: return self.pow(self.inv(a), -n)
        r = self.one()
        while n:
            if n & 1: r = self.mul(r, a)
            a = self.mul(a, a); n >>= 1
        return r
    def inv(self, a):
        # solve a * x = 1 via linear algebra over Q
        d = self.d
        M = [[self.mul(a, self.basis_vec(j))[i] for j in range(d)] for i in range(d)]
        # solve M x = e0
        A = [row[:] + [Fr(1) if i == 0 else Fr(0)] for i, row in enumerate(M)]
        for col in range(d):
            p = next(r for r in range(col, d) if A[r][col])
            A[col], A[p] = A[p], A[col]
            q = A[col][col]; A[col] = [x / q for x in A[col]]
            for r in range(d):
                if r != col and A[r][col]:
                    A[r] = [x - y * A[r][col] for x, y in zip(A[r], A[col])]
        return tuple(A[i][d] for i in range(d))
    def basis_vec(self, j):
        return tuple(Fr(1) if i == j else Fr(0) for i in range(self.d))
    def norm(self, a):
        # determinant of multiplication matrix
        d = self.d
        M = [[self.mul(a, self.basis_vec(j))[i] for j in range(d)] for i in range(d)]
        # det by Bareiss
        M = [row[:] for row in M]
        det = Fr(1)
        for col in range(d):
            p = next((r for r in range(col, d) if M[r][col]), None)
            if p is None: return Fr(0)
            if p != col:
                M[col], M[p] = M[p], M[col]; det = -det
            det *= M[col][col]
            for r in range(col + 1, d):
                if M[r][col]:
                    f = M[r][col] / M[col][col]
                    M[r] = [x - f * y for x, y in zip(M[r], M[col])]
        return det
    def eq(self, a, b): return all(x == y for x, y in zip(a, b))
    def eval(self, a, theta_val):
        t = mp.mpf(1); s = mp.mpf(0)
        for c in a:
            s += mp.mpf(c.numerator) / mp.mpf(c.denominator) * t
            t *= theta_val
        return s
    def roots(self, dps=100):
        if getattr(self, "_roots_cache", None) is not None and self._roots_dps >= dps:
            return self._roots_cache
        mp.mp.dps = max(dps, 60)
        coeffs = [mp.mpf(c.numerator) / mp.mpf(c.denominator) for c in reversed(self.m)]
        rts = mp.polyroots(coeffs, maxsteps=200, error=False)
        rts = [mp.re(r) for r in rts]
        rts.sort()
        self._roots_cache = rts
        self._roots_dps = dps
        return rts

def primes_above(m, p):
    """number of distinct primes above p in Q(theta) (p not dividing index), and factor degrees."""
    import sympy as sp
    x = sp.symbols('x')
    poly = sp.Poly([int(c) for c in reversed(m)], x, modulus=p)
    fac = poly.factor_list()
    # fac: (lc, [(f, mult), ...])
    degs = [(sp.degree(f[0]), f[1]) for f in fac[1]]
    return sum(mult for deg, mult in degs), degs
