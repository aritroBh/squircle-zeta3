"""sym2.py -- exact function symbols of GPLs G(a_1..a_n; 1) whose letters are rational functions of
parameters (sympy), with the Lambda^2 B_2 obstruction delta_{2,2}.  Letters of the symbol are the
irreducible factors (over Q) of the letter differences; constants are dropped.

Symbol recursion (Goncharov; as in A3's w4_symbol.py):
   S(G(x_1..x_n; e)) = sum_i S(G(..x_i omitted..; e)) (x) ( [x_i - x_{i-1}] - [x_i - x_{i+1}] ),
   x_0 = e (= 1), x_{n+1} = 0; coincident points contribute nothing (unit tangents).
"""
from collections import defaultdict
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp


class SymbolRing:
    def __init__(self, params):
        self.params = params
        self.alphabet = []          # list of irreducible polynomials (sympy expr, monic-ish)
        self.cache = {}

    def _key(self, f):
        return sp.srepr(sp.expand(f))

    def letter_index(self, f):
        f = sp.expand(f)
        for i, g in enumerate(self.alphabet):
            q = sp.cancel(f / g)
            if q.free_symbols.isdisjoint(self.params):
                return i
        self.alphabet.append(f)
        return len(self.alphabet) - 1

    def logs(self, expr):
        """dict letter_index -> exponent for a rational function (constants dropped)."""
        expr = sp.cancel(sp.together(expr))
        if expr == 0:
            return {}
        k = self._key(expr)
        if k in self.cache:
            return self.cache[k]
        num, den = sp.fraction(expr)
        out = defaultdict(int)
        for poly, sgn in ((num, 1), (den, -1)):
            _, facs = sp.factor_list(sp.expand(poly), *self.params)
            for fac, pw in facs:
                if fac.free_symbols.isdisjoint(self.params):
                    continue
                out[self.letter_index(fac)] += sgn * pw
        res = {i: v for i, v in out.items() if v}
        self.cache[k] = res
        return res

    def gpl_symbol(self, word, end=sp.Integer(1)):
        word = tuple(sp.sympify(w) for w in word)
        return self._gpl(word, sp.sympify(end))

    def _gpl(self, word, end):
        key = (tuple(sp.srepr(w) for w in word), sp.srepr(end))
        memo = self.__dict__.setdefault('_memo', {})
        if key in memo:
            return memo[key]
        if not word:
            res = {(): Fr(1)}
        else:
            pts = (end,) + word + (sp.Integer(0),)
            out = defaultdict(Fr)
            for j in range(1, len(word) + 1):
                x = pts[j]
                ell = defaultdict(Fr)
                d1 = sp.cancel(x - pts[j - 1]); d2 = sp.cancel(x - pts[j + 1])
                if d1 != 0:
                    for i, e in self.logs(d1).items():
                        ell[i] += e
                if d2 != 0:
                    for i, e in self.logs(d2).items():
                        ell[i] -= e
                ell = {i: v for i, v in ell.items() if v}
                if not ell:
                    continue
                sub = self._gpl(word[:j - 1] + word[j:], end)
                for w, c in sub.items():
                    for i, v in ell.items():
                        out[w + (i,)] += c * v
            res = {k: v for k, v in out.items() if v}
        memo[key] = res
        return res


def add(*ds, coeffs=None):
    out = defaultdict(Fr)
    for i, d in enumerate(ds):
        cf = Fr(1) if coeffs is None else Fr(coeffs[i])
        for k, v in d.items():
            out[k] += cf * v
    return {k: v for k, v in out.items() if v}


def _wedge(a, b):
    if a == b:
        return None, 0
    return ((a, b), 1) if a < b else ((b, a), -1)


def delta22(S):
    out = defaultdict(Fr)
    for word, v in S.items():
        if len(word) != 4:
            continue
        a, b, c, d = word
        p, s1 = _wedge(a, b)
        q, s2 = _wedge(c, d)
        if s1 and s2:
            key, s3 = _wedge(p, q)
            if s3:
                out[key] += v * s1 * s2 * s3
    return {k: v for k, v in out.items() if v}


def b2(ring, f):
    """Lambda^2 image of {f}_2 : (1-f) ^ f  as dict (i,j)->coeff, i<j."""
    out = defaultdict(Fr)
    A = ring.logs(1 - f); B = ring.logs(f)
    for i, a in A.items():
        for j, b in B.items():
            p, s = _wedge(i, j)
            if s:
                out[p] += a * b * s
    return {k: v for k, v in out.items() if v}


def wedge2(P, Q):
    """(sum P) ^ (sum Q) in Lambda^2(Lambda^2)"""
    out = defaultdict(Fr)
    for p, a in P.items():
        for q, b in Q.items():
            k, s = _wedge(p, q)
            if s:
                out[k] += a * b * s
    return {k: v for k, v in out.items() if v}
