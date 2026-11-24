"""Exact arithmetic in Q(zeta_N), zeta_N = e^{2 pi i / N}, power basis 1..zeta^{phi(N)-1},
elements as tuples of Fractions. Used for the exact Mobius-map enumeration (COR13C)."""
from fractions import Fraction as Fr
from functools import reduce
import sympy as sp

class Cyclo:
    def __init__(self, N):
        self.N = N
        x = sp.symbols('x')
        phi = sp.Poly(sp.cyclotomic_poly(N, x))
        self.d = d = phi.degree()
        self.mod = [Fr(int(c)) for c in phi.all_coeffs()]  # descending, monic, length d+1
        assert self.mod[0] == 1
        self.zero = (Fr(0),)*d
        self.one = (Fr(1),) + (Fr(0),)*(d-1)
        # roots zeta^k reduced, k = 0..N-1
        self.roots = []
        r = self.one
        zeta = (Fr(0), Fr(1)) + (Fr(0),)*(d-2) if d > 1 else self.one
        for k in range(N):
            self.roots.append(r)
            r = self.mul(r, zeta)
        assert r == self.one
    def add(self, x, y): return tuple(a+b for a,b in zip(x,y))
    def sub(self, x, y): return tuple(a-b for a,b in zip(x,y))
    def scal(self, x, s): return tuple(a*s for a in x)
    def mul(self, x, y):
        d = self.d
        tmp = [Fr(0)]*(2*d-1)
        for i,a in enumerate(x):
            if a:
                for j,b in enumerate(y):
                    if b: tmp[i+j] += a*b
        # reduce mod phi (monic, descending coeffs self.mod)
        for k in range(2*d-2, d-1, -1):
            c = tmp[k]
            if c:
                # subtract c * x^{k-d} * phi
                for i in range(1, d+1):
                    tmp[k-i] -= c*self.mod[i]
                tmp[k] = Fr(0)
        return tuple(tmp[:d])
    def is_root(self, x):
        """m such that x = zeta^m, else None"""
        for m, r in enumerate(self.roots):
            if x == r: return m
        return None

def mobius_images(C, i, j, k):
    """images of (0, 1/3, 1, -1, inf) under f with f(0)=z^i, f(1)=z^j, f(inf)=z^k.
    Returns tuple of ('root',m) / ('zero',) / ('inf',) for (1/3, -1)."""
    a, b, c = C.roots[i], C.roots[j], C.roots[k]
    bma = C.sub(b, a); cmb = C.sub(c, b)
    t1 = C.mul(c, bma); t2 = C.mul(a, cmb)
    for num, den in ((C.add(t1, C.scal(t2, Fr(3))), C.add(bma, C.scal(cmb, Fr(3)))),
                     (C.sub(t2, t1), C.sub(cmb, bma))):
        if den == C.zero:
            yield ('inf',)
        elif num == C.zero:
            yield ('zero',)
        else:
            hit = None
            for m, r in enumerate(C.roots):
                if C.mul(r, den) == num:
                    hit = ('root', m); break
            yield hit if hit else ('other',)
