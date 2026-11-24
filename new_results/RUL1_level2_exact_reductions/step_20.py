import subprocess
code = r'''
"""rul1.py -- exact level-2 MZV (Euler sum) regularized-double-shuffle engine, weights <= 5.

Words: tuples over alphabet {0, 1, -1}, letter a = form dt/(t-a); value G(w;1) (iterated
integral, 0<t_1<...<t_n<1, first letter outermost). Convergent iff first != 1 and last != 0.

Stuples: tuples of (s, eps), eps in {+1,-1}; value zeta(S) = sum_{n_1>...>n_d>=1} prod eps_i^{n_i}/n_i^{s_i}.
Theorem (standard conversion): zeta(S) = (-1)^d G(w(S);1), where w(S) has blocks 0^{s_i-1} e_i,
e_i = eps_1*...*eps_i. All relations used:
  DS   : finite double shuffle  u ш v (words) == st(u) * st(v) (stuples, transported)   [proved: absolute convergence]
  RDS  : shuffle-reg value of word(S) == stuffle-reg value of S (divergent S)            [proved: Ihara-Kaneko-Zagier comparison at T=0]
  DIST : distribution  sum_{delta in {±1}^d} zeta(s_i, delta_i eps_i) = 2^{d-w} zeta(s_i, +1), only for s_1 >= 2 [proved: absolute convergence]
  DUAL : duality t -> (1-t)/(1+t) on convergent words                                     [proved: change of variables]
"""
from fractions import Fraction as Fr
from functools import cache
import math, itertools

# ---------------- basic algebra ----------------

def stuple_to_word(S):
    w = []
    e = 1
    for (s, eps) in S:
        e *= eps
        w += [0]*(s-1) + [e]
    return tuple(w)

def word_to_stuple(w):
    S = []
    e_prev = 1
    run = 0
    for a in w:
        if a == 0:
            run += 1
        else:
            S.append((run+1, a*e_prev))
            e_prev = a
            run = 0
    if run:
        raise ValueError("word ends in 0")
    return tuple(S)

def depth_word(w):
    return sum(1 for a in w if a != 0)

def convergent_word(w):
    return (not w) or (w[0] != 1 and w[-1] != 0)

def convergent_stuple(S):
    return (not S) or S[0] != (1, 1)

@cache
def shuffle(u, v):
    if not u: return {v: Fr(1)}
    if not v: return {u: Fr(1)}
    out = {}
    for w, c in shuffle(u[1:], v).items():
        key = (u[0],)+w
        out[key] = out.get(key, Fr(0)) + c
    for w, c in shuffle(u, v[1:]).items():
        key = (v[0],)+w
        out[key] = out.get(key, Fr(0)) + c
    return out

@cache
def stuffle(S, T):
    if not S: return {T: Fr(1)}
    if not T: return {S: Fr(1)}
    out = {}
    for w, c in stuffle(S[1:], T).items():
        key = (S[0],)+w
        out[key] = out.get(key, Fr(0)) + c
    for w, c in stuffle(S, T[1:]).items():
        key = (T[0],)+w
        out[key] = out.get(key, Fr(0)) + c
    m = ((S[0][0]+T[0][0], S[0][1]*T[0][1]),)
    for w, c in stuffle(S[1:], T[1:]).items():
        key = m + w
        out[key] = out.get(key, Fr(0)) + c
    return out

def sh_mul(A, B):
    """product of two word-combos via shuffle (G-value product)."""
    out = {}
    for u, a in A.items():
        for v, b in B.items():
            for w, c in shuffle(u, v).items():
                out[w] = out.get(w, Fr(0)) + a*b*c
    return {k: v for k, v in out.items() if v}

def st_mul(A, B):
    out = {}
    for u, a in A.items():
        for v, b in B.items():
            for w, c in stuffle(u, v).items():
                out[w] = out.get(w, Fr(0)) + a*b*c
    return {k: v for k, v in out.items() if v}

# ---------------- regularizations ----------------

@cache
def _strip_lead1(w):
    """w = sum c * (1^a ш u), u not starting with 1.  Returns dict (a,u)->c."""
    w = tuple(w)
    if not w or w[0] != 1: return {(0, w): Fr(1)}
    t = 0
    while t < len(w) and w[t] == 1: t += 1
    v = w[t:]
    if not v: return {(t, ()): Fr(1, math.factorial(t))}
    w2 = (1,)*(t-1) + v
    P = dict(shuffle((1,), w2))
    cw = P.pop(w)
    assert cw == t, (w, cw, t)
    out = {}
    for (a, u), c in _strip_lead1(w2).items():
        key = (a+1, u)
        out[key] = out.get(key, Fr(0)) + Fr(1, t)*c*(a+1)
    for T, cT in P.items():
        for (a, u), c in _strip_lead1(T).items():
            key = (a, u)
            out[key] = out.get(key, Fr(0)) - Fr(1, t)*cT*c
    return {k: v for k, v in out.items() if v}

@cache
def _strip_trail0(w):
    """w = sum c * (u ш 0^b), u not ending with 0.  Returns dict (b,u)->c."""
    w = tuple(w)
    if not w or w[-1] != 0: return {(0, w): Fr(1)}
    t = 0
    while t < len(w) and w[-1-t] == 0: t += 1
    v = w[:-t]
    if not v: return {(t, ()): Fr(1, math.factorial(t))}
    w2 = v + (0,)*(t-1)
    P = dict(shuffle(w2, (0,)))
    cw = P.pop(w)
    assert cw == t, (w, cw, t)
    out = {}
    for (b, u), c in _strip_trail0(w2).items():
        key = (b+1, u)
        out[key] = out.get(key, Fr(0)) + Fr(1, t)*c*(b+1)
    for T, cT in P.items():
        for (b, u), c in _strip_trail0(T).items():
            key = (b, u)
            out[key] = out.get(key, Fr(0)) - Fr(1, t)*cT*c
    return {k: v for k, v in out.items() if v}

@cache
def sh_reg(w):
    """shuffle-regularized G-value (T=0) as dict convergent-word -> Fraction."""
    w = tuple(w)
    out = {}
    for (a, u), c in _strip_lead1(w).items():
        if a != 0: continue
        for (b, v), c2 in _strip_trail0(u).items():
            if b != 0: continue
            out[v] = out.get(v, Fr(0)) + c*c2
    return {k: v for k, v in out.items() if v}

@cache
def st_reg(S):
    """stuffle-regularized zeta-value (T=0) as dict convergent-stuple -> Fraction."""
    S = tuple(S)
    if not S or S[0] != (1, 1): return {S: Fr(1)}
    Sp = S[1:]
    P = dict(stuffle(((1, 1),), Sp))
    cS = P.pop(S)
    assert cS == 1
    out = {}
    for T, cT in P.items():
        for U, cU in st_reg(T).items():
            out[U] = out.get(U, Fr(0)) - cT*cU
    return {k: v for k, v in out.items() if v}

def stuple_set_to_words(d):
    """transport a stuple-combo to a convergent-word combo via zeta(S)=(-1)^d G(w(S))."""
    out = {}
    for S, c in d.items():
        w = stuple_to_word(S)
        cw = sh_reg(w) if not convergent_word(w) else {w: Fr(1)}
        for v, c2 in cw.items():
            out[v] = out.get(v, Fr(0)) + c*c2
    return {k: v for k, v in out.items() if v}

# note: stuple_set_to_words uses sh_reg for divergent words; for convergent stuples the word is
# convergent already. The sign (-1)^d is handled by the conversion theorem applied uniformly:
# we write the VALUE map as G on words; zeta(S) = (-1)^d G(w(S)). In DS relations the signs
# cancel (uniform depth), see relations() below where we transport stuffle products to words.

# ---------------- word/stuple enumeration ----------------

def words_of_weight(w):
    return list(itertools.product((0, 1, -1), repeat=w)) if w > 0 else [()]

def conv_words_of_weight(w):
    return [u for u in words_of_weight(w) if convergent_word(u)]

def stuples_of_weight(w):
    out = []
    def rec(rem, acc):
        if rem == 0:
            out.append(tuple(acc)); return
        for s in range(1, rem+1):
            for e in (1, -1):
                rec(rem-s, acc+[(s, e)])
    rec(w, [])
    return out

def conv_stuples_of_weight(w):
    return [S for S in stuples_of_weight(w) if convergent_stuple(S)]

# ---------------- relations ----------------
# A relation is dict convergent-word -> Fraction (meaning: this combination has value 0).

def relation_DS(u, v):
    """finite double shuffle for convergent u,v: u ш v - st(u)*st(v) transported to words."""
    lhs = shuffle(u, v)
    Su, Sv = word_to_stuple(u), word_to_stuple(v)
    rhs = stuple_set_to_words(stuffle(Su, Sv))
    out = dict(lhs)
    for w, c in rhs.items():
        out[w] = out.get(w, Fr(0)) - c
    return {k: x for k, x in out.items() if x}

def relation_RDS(S):
    """regularized double shuffle for divergent stuple S: sh_reg(w(S)) - st_reg(S) transported."""
    w = stuple_to_word(S)
    lhs = sh_reg(w)
    rhs = stuple_set_to_words(st_reg(S))
    out = dict(lhs)
    for ww, c in rhs.items():
        out[ww] = out.get(ww, Fr(0)) - c
    return {k: x for k, x in out.items() if x}

def relation_DIST(S):
    """distribution for convergent S with s_1 >= 2: sum_{delta} zeta(s_i, delta_i eps_i) = 2^{d-w} zeta(s_i,+1)."""
    d = len(S)
    w = sum(s for s, e in S)
    out = {}
    for deltas in itertools.product((1, -1), repeat=d):
        T = tuple((s, dd*e) for (s, e), dd in zip(S, deltas))
        out = _addinto(out, stuple_set_to_words(st_reg(T) if not convergent_stuple(T) else {T: Fr(1)}), Fr(1))
    plain = tuple((s, 1) for s, e in S)
    out = _addinto(out, stuple_set_to_words(st_reg(plain) if not convergent_stuple(plain) else {plain: Fr(1)}), -Fr(2**(d-w)))
    return {k: x for k, x in out.items() if x}

def _addinto(A, B, c):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, Fr(0)) + v*c
    return {k: x for k, x in out.items() if x}

def dual_word(w):
    """t -> (1-t)/(1+t): letters pull back as 0 -> (1)-(-1), 1 -> (0)-(-1), -1 -> -(-1); reverse word."""
    letter = {0: ((1, Fr(1)), (-1, Fr(-1))),
              1: ((0, Fr(1)), (-1, Fr(-1))),
              -1: ((-1, Fr(-1)),)}
    out = {(): Fr(1)}
    for a in reversed(w):
        new = {}
        for suff, c in out.items():
            for l, sgn in letter[a]:
                nw = (l,) + suff
                new[nw] = new.get(nw, Fr(0)) + c*sgn
        out = new
    return {k: v for k, v in out.items() if v}

def relation_DUAL(w):
    """duality for convergent w: G(w) = G(dual(w))."""
    out = {w: Fr(1)}
    for v, c in dual_word(w).items():
        out[v] = out.get(v, Fr(0)) - c
    return {k: x for k, x in out.items() if x}

def relations_at_weight(w, families=('DS', 'RDS', 'DIST', 'DUAL')):
    rels = []
    if 'DS' in families:
        cw = [u for u in conv_words_of_weight_less(w)]
        for i, u in enumerate(cw):
            for v in cw[i:]:
                if len(u) + len(v) != w: continue
                if len(u) == 0 or len(v) == 0: continue
                rels.append(('DS', (u, v), relation_DS(u, v)))
    if 'RDS' in families:
        for S in stuples_of_weight(w):
            if not convergent_stuple(S):
                rels.append(('RDS', S, relation_RDS(S)))
    if 'DIST' in families:
        for S in conv_stuples_of_weight(w):
            if S and S[0][0] >= 2:
                rels.append(('DIST', S, relation_DIST(S)))
    if 'DUAL' in families:
        for u in conv_words_of_weight(w):
            rels.append(('DUAL', u, relation_DUAL(u)))
    return rels

@cache
def conv_words_of_weight_less(w):
    out = []
    for k in range(1, w):
        out += conv_words_of_weight(k)
    return tuple(out)
'''
import pathlib
pathlib.Path('rul1.py').write_text(code)
print("written", len(code.splitlines()), "lines")
