
# step_13.py -- symbol via Goncharov's coproduct formula (third method); test at c0=1/2
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
from functools import lru_cache
import sympy as sp, json, re

EPS = sp.Symbol('EPSREG3')
def cn(e):
    e = sp.simplify(e)
    if e is EPS: return e
    ne = sp.simplify(-e)
    ke, kn = (str(e).startswith('-'), str(e)), (str(ne).startswith('-'), str(ne))
    return e if ke <= kn else ne
def isz(e):
    e = sp.simplify(e)
    return e == 0 or e in (1, -1, sp.I, -sp.I)
def dl(P, Q):
    d = sp.simplify(P[1] - Q[1])
    if d == 0:
        if P[0] == 'c' and Q[0] == 'c': raise AssertionError("collide")
        return EPS
    return cn(d)

@lru_cache(maxsize=None)
def symC(a0, letters, aN):
    """Symbol of regularized I(a0; a_1..a_n; aN) via Goncharov's coproduct recursion:
       S = sum_i ([aN - a_i] - [a0 - a_i]) (x) S( I(a0; a_1..a_{i-1}; a_i) shuffle I(a_i; a_{i+1}..a_n; aN) ).
       a0, aN: ('c', expr); letters perturbed a_i + eps^i."""
    n = len(letters)
    A0 = ('c', a0); AN = ('c', aN)
    LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
    if n == 0:
        return {(): 1}
    out = {}
    for i in range(n):
        ai = LL[i]
        e1 = dl(AN, ai)   # [aN - a_i]
        e0 = dl(A0, ai)   # [a0 - a_i]
        SL = symC(a0, tuple(l[1] for l in LL[:i]), ai[1]) if i > 0 else {(): 1}
        SR = symC(ai[1], tuple(l[1] for l in LL[i+1:]), aN) if i < n - 1 else {(): 1}
        Sseg = sym_shuffle(SL, SR)
        for sgn, ent in ((+1, e1), (-1, e0)):
            if ent is EPS or isz(ent): continue
            for w2, c2 in Sseg.items():
                out[(ent,) + w2] = out.get((ent,) + w2, 0) + sgn * c2
    return {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}

x = sp.Symbol('x')
print("C: S(G(0,1;x)) =", symC(sp.Integer(0), (0, 1), x))
print("mine:", sym_G((0,1), x))
print("C: S(G(1,0;x)) =", symC(sp.Integer(0), (1, 0), x))
print("mine:", sym_G((1,0), x))
print("C: S(G(0,0,0,1;x)) =", symC(sp.Integer(0), (0,0,0,1), x))
print("mine:", sym_G((0,0,0,1), x))
print("C: S(Li_2) shuffle tests:")
def check_shuffleC(u, v):
    lhs = sym_shuffle(symC(sp.Integer(0), u, x), symC(sp.Integer(0), v, x))
    rhs = {}
    for w in shuffle_words(u, v):
        rhs = sym_add(rhs, symC(sp.Integer(0), w, x))
    return sym_sub(lhs, rhs) == {}
pairs = [((0,),(0,)), ((1,),(-1,)), ((0,),(1,)), ((0,1),(0,1)), ((0,),(0,1)), ((1,-1),(0,1)),
         ((0,1),(0,-1)), ((0,0,1),(0,)), ((0,0,1),(1,)), ((1,),(0,0,1)), ((0,),(0,0,0,1)),
         ((0,1),(1,-1)), ((1,),(1,0,-1)), ((-1,),(0,1,-1)), ((0,-1),(0,0,1))]
ok = all(check_shuffleC(u, v) for u, v in pairs)
print("  coproduct-recursion shuffle checks:", "PASS" if ok else "FAIL")
print("C: S(G(0,0,-1,1;x)) =", symC(sp.Integer(0), (0,0,-1,1), x))
print("mine (same word)   =", sym_G((0,0,-1,1), x))
