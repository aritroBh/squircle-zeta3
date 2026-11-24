"""Own implementation of formal MZV-word algebra (alphabet {0,1,-1}), independent of THX1 code.
G(a1..an;1) = int_{1>t1>..>tn>0} prod dt_i/(t_i-a_i).  Convergent iff a1 != 1 and an != 0.
Sum dictionary: zeta(s1..sd; e1..ed) = sum_{n1>..>nd>=1} prod e_i^{n_i}/n_i^{s_i} = (-1)^d G(0^{s1-1} b1 .. 0^{sd-1} bd), b_j = e_1...e_j.
"""
from functools import lru_cache
from fractions import Fraction as F

def add(acc, d, c=1):
    for k, v in d.items():
        acc[k] = acc.get(k, 0) + c*v

@lru_cache(None)
def _sh(u, v):
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = {}
    for w, c in _sh(u[1:], v).items():
        k = (u[0],)+w; out[k] = out.get(k, 0)+c
    for w, c in _sh(u, v[1:]).items():
        k = (v[0],)+w; out[k] = out.get(k, 0)+c
    return out

def shuffle(x, y):
    """bilinear shuffle product of dicts word->coef"""
    out = {}
    for u, a in x.items():
        for v, b in y.items():
            for w, c in _sh(u, v).items():
                out[w] = out.get(w, 0) + a*b*c
    return {k: v for k, v in out.items() if v}

def word_to_sum(w):
    """convergent-at-0 word (last letter !=0) -> tuple of (s,e)"""
    assert w and w[-1] != 0
    S = []; z = 0; prev = 1
    for a in w:
        if a == 0: z += 1
        else:
            S.append((z+1, a*prev)); prev = a; z = 0
    return tuple(S)

def sum_to_word(S):
    w = []; b = 1
    for s, e in S:
        b *= e; w += [0]*(s-1)+[b]
    return tuple(w)

@lru_cache(None)
def _st(a, b):
    if not a: return {b: 1}
    if not b: return {a: 1}
    out = {}
    for t, c in _st(a[1:], b).items():
        k = (a[0],)+t; out[k] = out.get(k, 0)+c
    for t, c in _st(a, b[1:]).items():
        k = (b[0],)+t; out[k] = out.get(k, 0)+c
    m = (a[0][0]+b[0][0], a[0][1]*b[0][1])
    for t, c in _st(a[1:], b[1:]).items():
        k = (m,)+t; out[k] = out.get(k, 0)+c
    return out

def stuffle_words(u, v):
    """G(u)G(v) via stuffle of sums: returns dict word->coef"""
    a, b = word_to_sum(u), word_to_sum(v)
    sg = (-1)**(len(a)+len(b))
    out = {}
    for t, c in _st(a, b).items():
        w = sum_to_word(t)
        out[w] = out.get(w, 0) + sg*(-1)**len(t)*c
    return out

def DS(u, v):
    """finite double shuffle relation sh(u,v)-st(u,v) for convergent words u,v (vector word->int)"""
    out = {}
    add(out, _sh(u, v)); add(out, stuffle_words(u, v), -1)
    return {k: x for k, x in out.items() if x}

_g = {0: ((-1, 1), (1, -1)), 1: ((-1, 1), (0, -1)), -1: ((-1, 1),)}
def dual_image(w):
    terms = {(): 1}
    out = {}
    def rec(i, acc, c):
        if i == len(w):
            k = tuple(reversed(acc)); out[k] = out.get(k, 0)+c; return
        for b, cf in _g[w[i]]:
            rec(i+1, acc+[b], c*cf)
    rec(0, [], 1)
    return out
def DUAL(w):
    out = {w: 1}; add(out, dual_image(w), -1)
    return {k: x for k, x in out.items() if x}

def conv(w): return bool(w) and w[0] != 1 and w[-1] != 0

def zeta_word(n): return {(0,)*(n-1)+(1,): F(-1)}
def zbar(S):  # zeta of sum tuple S -> word dict, with (-1)^d
    return {sum_to_word(S): F((-1)**len(S))}
def mono(*fs):
    out = {(): F(1)}
    for f in fs: out = shuffle(out, f)
    return out
