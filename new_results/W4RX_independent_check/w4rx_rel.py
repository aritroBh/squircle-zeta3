"""W4RX relation families, implemented from the definitions (imports only w4rx_core).

All computations are in G notation (see w4rx_core).  Enumeration ORDER and the row
normalisation follow the notes' specification (stated in result.md) so that the certificate's
row labels can be read; the CONTENT of every row is computed here from the mathematics.
"""
from w4rx_core import *
from functools import lru_cache
from itertools import product
from math import factorial

# ---------------------------------------------------------------- shuffle regularisation
@lru_cache(maxsize=None)
def reg(w):
    """Shuffle regularisation at the upper end point: the unique shuffle-algebra map with
    reg(word)=word for words not starting with letter 0 (the point 1) and reg((0,))=0.
    For w = 1^m w' (m>=1 leading letters 0), (0) sh 1^(m-1)w' = m*w + (words with m-1 leading
    zeros), and reg((0) sh x) = reg((0)) reg(x) = 0."""
    if not w or w[0] != 0: return {w: F(1)}
    m = 0
    while m < len(w) and w[m] == 0: m += 1
    terms = dict(sh((0,), w[1:]))
    c = terms.pop(w)
    assert c == m
    out = {}
    for t, a in terms.items():
        out = vadd((out, 1), (reg(t), F(-a, c)))
    return out
def regv(vec):
    out = {}
    for w, c in vec.items(): out = vadd((out, 1), (reg(w), c))
    return out

# ---------------------------------------------------------------- series <-> G words
def series_of(w):
    """G(0^(s1-1),c1,...,0^(sd-1),cd;1) = (-1)^d Li_{s1..sd}(1/c1, c1/c2, ..., c(d-1)/cd).
    Roots of unity are stored by exponents: x_i = zeta^(e_i)."""
    out = []; s = 1; prev = 0
    for a in w:
        if a == ZERO: s += 1
        else:
            out.append((s, (prev - a) % 8)); prev = a; s = 1
    assert s == 1, 'trailing letter 0'
    return tuple(out)
def word_of(series):
    out = []; cum = 0
    for s, e in series:
        cum = (cum + e) % 8                  # b_i = x1...xi = zeta^cum, letter c_i = 1/b_i
        out.extend([ZERO] * (s - 1) + [(-cum) % 8])
    return tuple(out)
@lru_cache(maxsize=None)
def stuffle(u, v):
    """Product of the series sum_{n1>..} x^n/n^s: quasi-shuffle with (s,x)(t,y) -> (s+t, xy)."""
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = Counter()
    for w, c in stuffle(u[1:], v).items(): out[(u[0],) + w] += c
    for w, c in stuffle(u, v[1:]).items(): out[(v[0],) + w] += c
    m = (u[0][0] + v[0][0], (u[0][1] + v[0][1]) % 8)
    for w, c in stuffle(u[1:], v[1:]).items(): out[(m,) + w] += c
    return dict(out)

def double_shuffle(u, v):
    """reg( stuffle(u,v) - shuffle(u,v) ), u any word without trailing 0, v admissible;
    the stuffle side in G-words: G(u)G(v) = (-1)^(du+dv) Li(u)Li(v) = sum_t (-1)^(du+dv+dt) G(t)."""
    su, sv = series_of(u), series_of(v)
    st = {}
    for t, c in stuffle(su, sv).items():
        w = word_of(t)
        st[w] = st.get(w, 0) + c * (-1) ** (len(su) + len(sv) + len(t))
    return regv(vadd((st, 1), (sh(u, v), -1)))

def distribution(w, d):
    """w: G-word with letters in {0 (ZERO)} U {zeta^m : d | m}, last letter nonzero.
    t = s^d:  G(w;1) = d^{#0} sum_lifts G(lift;1); regularised:
    d^{#0} sum_lifts reg(lift) = sum_{k<=m} C^{sh k}/k! sh reg(w[k:]),  m = number of leading
    letters 1 (letter 0), C = sum_{r^d=1, r!=1} G(r;1) = log d."""
    C = {((8 // d) * k,): F(1) for k in range(1, d)}
    lifts = [[(ZERO,)] if a == ZERO else [((a + 8 * k) // d % 8,) for k in range(d)] for a in w]
    left = Counter()
    for tup in product(*lifts): left[tuple(x[0] for x in tup)] += d ** w.count(ZERO)
    right = dict(reg(w)); Ck = {(): F(1)}
    m = 0
    while m < len(w) and w[m] == 0: m += 1
    for k in range(1, m + 1):
        Ck = vmul(Ck, C)
        right = vadd((right, 1), (vmul(Ck, reg(w[k:])), F(1, factorial(k))))
    return regv(vadd((dict(left), 1), (right, -1)))

def branch_w1(m):
    """G(zeta^m;1)-G(zeta^-m;1) = log(-zeta^-m) = i pi (1-m/4), and G(i)-G(-i) = i pi/2:
    2(G(zeta^m)-G(zeta^-m)) - (4-m)(G(zeta^2)-G(zeta^6)) = 0."""
    return vadd(({(m,): 2}, 1), ({((-m) % 8,): 2}, -1), ({(2,): 1}, -(4 - m)), ({(6,): 1}, 4 - m))

# ---------------------------------------------------------------- enumeration (notes' order)
def h_words(n, admissible_only=True):
    if not n: return ((),)
    return tuple(w for w in product(range(-1, 8), repeat=n)
                 if w[-1] != -1 and (not admissible_only or w[0] != 0))

def standard_rows(n, lower, log=print):
    """Rows in the notes' enumeration order: weight-1 branch rows (n=1) or double shuffle
    (a=1..n-1, u in all H-words of weight a, v admissible of weight b), then distribution
    d=2,4,8, then lifts of lower bases; primitive normalisation; first occurrence kept."""
    rows = []; seen = {}; fam = []
    def emit(gvec, family):
        r = primitive_H(to_H(gvec))
        if r and r not in seen:
            seen[r] = len(rows); rows.append(r); fam.append(family)
    if n == 1:
        for j in range(1, 8):
            m = (-j) % 8
            emit(branch_w1(m), ('branch', j))
    else:
        for a in range(1, n):
            b = n - a
            for u in h_words(a, False):
                ua = u[0] != 0
                if ua and a > b: continue
                for v in h_words(b):
                    if ua and a == b and u > v: continue
                    emit(double_shuffle(g_of_h(u)[0], g_of_h(v)[0]), ('ds', u, v))
    for d in (2, 4, 8):
        for w in product((-1,) + tuple(range(0, 8, d)), repeat=n):
            if w[-1] == -1: continue
            emit(distribution(g_of_h(w)[0], d), ('dist', d, w))
    for k, basis in lower.items():
        for row in basis:
            for w in h_words(n - k):
                emit(vmul(from_H(dict(row)), {g_of_h(w)[0]: 1}), ('lift', k, w))
    return rows, fam
