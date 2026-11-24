"""hpl_tools.py -- exact toolkit for harmonic polylogarithms H_w(c) with letters {0, 1, -1}.

Conventions (Remiddi-Vermaseren, as in G1):  H_{()} = 1,  H_{a,w}(c) = int_0^c f_a(t) H_w(t) dt,
   f_0 = 1/t,  f_1 = 1/(1-t),  f_{-1} = 1/(1+t).
So H_1 = -log(1-c), H_{-1} = log(1+c), H_{0,0,0,1} = Li_4(c), H_{0,..,0,-1} = -Li_n(-c).
Linear combinations: dict word(tuple) -> Fraction.  Only words with nonzero LAST letter are used
(these are analytic at c = 0 and vanish there), plus the empty word.

Also: GPL -> HPL(c) conversion for the "r-family" of the squircle integral (see a04_rfamily.py):
GPLs G(a_1..a_n; 1) with letters among symbolic points whose positions are functions of c.
The total differential (Goncharov):
   d G(a_1..a_n; a_0) = sum_i G(..a_i omitted..; a_0) [ dlog(a_{i-1}-a_i) - dlog(a_{i+1}-a_i) ],
   a_0 = endpoint (=1), a_{n+1} = 0, terms with coincident points dropped (unit-tangent rule).
If every dlog(difference) is a Q-combination of dlog c, dlog(1-c), dlog(1+c) (+ constants) and every
GPL tends to 0 as c -> 0, then G = iterated integral from 0, i.e. a Q-combination of H_w(c).
"""
from fractions import Fraction as Fr
from collections import defaultdict
from functools import lru_cache
import itertools

LETTERS = (0, 1, -1)


def clean(d):
    return {k: v for k, v in d.items() if v != 0}


def add(*ds, coeffs=None):
    out = defaultdict(Fr)
    for i, d in enumerate(ds):
        cf = Fr(1) if coeffs is None else Fr(coeffs[i])
        for k, v in d.items():
            out[k] += cf * v
    return clean(out)


def scale(d, s):
    s = Fr(s)
    return clean({k: v * s for k, v in d.items()})


@lru_cache(maxsize=None)
def shuffle_words(u, v):
    """multiset of shuffles of tuples u, v as dict word->count"""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = defaultdict(int)
    for w, c in shuffle_words(u[1:], v).items():
        out[(u[0],) + w] += c
    for w, c in shuffle_words(u, v[1:]).items():
        out[(v[0],) + w] += c
    return dict(out)


def shuffle(A, B):
    out = defaultdict(Fr)
    for u, a in A.items():
        for v, b in B.items():
            for w, c in shuffle_words(u, v).items():
                out[w] += a * b * c
    return clean(out)


def prepend(a, A):
    return {(a,) + w: v for w, v in A.items()}


def weight_part(A, k):
    return {w: v for w, v in A.items() if len(w) == k}


# ---------------------------------------------------------------- Taylor series / numerics
def hpl_taylor(word, N, ring=Fr):
    """Taylor coefficients [c^0..c^N] of H_word(c) (exact Fractions by default)."""
    if len(word) == 0:
        s = [ring(0)] * (N + 1)
        s[0] = ring(1)
        return s
    inner = hpl_taylor(word[1:], N, ring)
    a = word[0]
    out = [ring(0)] * (N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m] / m
    else:
        acc = ring(0)
        for m in range(0, N):
            acc = acc * a + inner[m]
            out[m + 1] = acc / (m + 1)
    return out


def taylor_of(A, N, ring=Fr):
    out = [ring(0)] * (N + 1)
    cache = {}
    for w, v in A.items():
        if w not in cache:
            cache[w] = hpl_taylor(w, N, ring)
        s = cache[w]
        for i in range(N + 1):
            out[i] += ring(v) * s[i] if ring is not Fr else v * s[i]
    return out


def hpl_numeric_all(words, x, nterms):
    """Evaluate H_w(x) for all words (mpmath), via truncated Taylor series (|x| < 1).
    Uses memoized nested evaluation on suffixes."""
    import mpmath as mp
    memo = {}

    def ser(w):
        if w in memo:
            return memo[w]
        if len(w) == 0:
            s = [mp.mpf(0)] * (nterms + 1)
            s[0] = mp.mpf(1)
        else:
            inner = ser(w[1:])
            a = w[0]
            s = [mp.mpf(0)] * (nterms + 1)
            if a == 0:
                for m in range(1, nterms + 1):
                    s[m] = inner[m] / m
            else:
                acc = mp.mpf(0)
                for m in range(0, nterms):
                    acc = acc * a + inner[m]
                    s[m + 1] = acc / (m + 1)
        memo[w] = s
        return s
    vals = {}
    for w in words:
        s = ser(w)
        tot = mp.mpf(0)
        p = mp.mpf(1)
        for m in range(nterms + 1):
            tot += s[m] * p
            p *= x
        vals[w] = tot
    return vals


def eval_lin(A, x, nterms):
    import mpmath as mp
    vals = hpl_numeric_all(list(A.keys()), x, nterms)
    return mp.fsum(mp.mpf(v.numerator) / v.denominator * vals[w] for w, v in A.items())


# ---------------------------------------------------------------- symbols / obstruction
IDX = {0: 0, 1: 1, -1: 2}          # symbol letters: 0=[c], 1=[1-c], 2=[1+c]
SGN = {0: 1, 1: -1, -1: 1}


def symbol(A):
    """function symbol (constants dropped) in letters 0=[c],1=[1-c],2=[1+c];
    S(H_{a1..an}) = s(an) x ... x s(a1)."""
    out = defaultdict(Fr)
    for w, v in A.items():
        sg = 1
        for a in w:
            sg *= SGN[a]
        out[tuple(IDX[a] for a in reversed(w))] += v * sg
    return clean(out)


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
    return clean(out)


# ---------------------------------------------------------------- GPL -> HPL(c)
class GPLConverter:
    """Points: dict name -> None (symbolic).  dlogdiff(x, y): dict {0,1,-1} -> Fraction giving
    d/dc log(x - y) = sum_a coeff_a f_a(c)  (f_0 = 1/c, f_1 = 1/(1-c), f_-1 = 1/(1+c)),
    i.e. dlog c -> f_0,  dlog(1-c) -> -f_1,  dlog(1+c) -> f_-1.  Must be symmetric in x, y.
    Endpoint name 'E' (value 1), base point name 'B' (value 0)."""

    def __init__(self, dlogdiff):
        self.dlogdiff = dlogdiff
        self.memo = {}

    def D(self, x, y):
        if x == y:
            return {}
        return self.dlogdiff(x, y)

    def conv(self, word):
        """HPL(c)-expansion of G(word; 1); word: tuple of point names; G(word;1) -> 0 as c -> 0 assumed."""
        if word in self.memo:
            return self.memo[word]
        if len(word) == 0:
            res = {(): Fr(1)}
        else:
            n = len(word)
            pts = ('E',) + tuple(word) + ('B',)
            acc = defaultdict(Fr)
            for i in range(1, n + 1):
                form = defaultdict(Fr)
                for a, v in self.D(pts[i - 1], pts[i]).items():
                    form[a] += v
                for a, v in self.D(pts[i + 1], pts[i]).items():
                    form[a] -= v
                form = clean(form)
                if not form:
                    continue
                sub = self.conv(word[:i - 1] + word[i:])
                for a, fv in form.items():
                    for w, sv in sub.items():
                        acc[(a,) + w] += fv * sv
            res = clean(acc)
        self.memo[word] = res
        return res


def gpl_shuffle_power(terms, k):
    """terms: dict word(tuple of point names) -> Fraction (a linear combination of GPLs in x,
    including the empty word for constants);  returns the k-th shuffle power."""
    res = {(): Fr(1)}
    for _ in range(k):
        res = shuffle(res, terms)
    return res
