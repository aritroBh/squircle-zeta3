t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
from functools import lru_cache
import sympy as sp, json, re, itertools

EPS = sp.Symbol('EPSREG2')
def canon6(e):
    e = sp.simplify(e)
    if e is EPS: return e
    ne = sp.simplify(-e)
    ke, kn = (str(e).startswith('-'), str(e)), (str(ne).startswith('-'), str(ne))
    return e if ke <= kn else ne
def iszero6(e):
    e = sp.simplify(e)
    return e == 0 or e in (1, -1, sp.I, -sp.I)
def dlog6(P, Q):
    d = sp.simplify(P[1] - Q[1])
    if d == 0:
        if P[0] == 'c' and Q[0] == 'c': raise AssertionError("collide")
        return EPS
    return canon6(d)
@lru_cache(maxsize=None)
def sym6(a0, letters, aN):
    n = len(letters)
    A0 = ('c', a0); AN = ('c', aN)
    LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
    if n == 0: return {(): 1}
    out = {}
    for j in range(n):
        p = A0 if j == 0 else LL[j - 1]
        q = AN if j == n - 1 else LL[j + 1]
        num = dlog6(q, LL[j]); den = dlog6(p, LL[j])
        Ssub = sym6(a0, letters[:j] + letters[j+1:], aN)
        for sgn, ent in ((+1, num), (-1, den)):
            if ent is EPS or iszero6(ent): continue
            for w2, c2 in Ssub.items():
                out[(ent,) + w2] = out.get((ent,) + w2, 0) + sgn * c2
    return {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}

def decomposable(T, primes):
    basis = list(itertools.product(primes, repeat=4))
    idx = {m: i for i, m in enumerate(basis)}
    rows = []
    for u in (1, 2, 3):
        for m1 in itertools.product(primes, repeat=u):
            for m2 in itertools.product(primes, repeat=4 - u):
                r = [0]*len(basis)
                for w in shuffle_words(m1, m2): r[idx[w]] += 1
                rows.append(r)
    P = sp.Matrix(rows)
    dv = sp.Matrix([T.get(m, 0) for m in basis])
    return P.T.rank() == P.T.row_join(dv).rank(), P.rank()

x = sp.Symbol('x')
primes_xy = None
# check per-word difference decomposability at 1/2 (primes {2,3})
c0 = sp.Rational(1, 2)
words = [(0,0,-1,1), (0,0,1,-1), (1,-1,0,1), (-1,1,0,-1), (1,0,-1,1), (-1,1,1,-1)]
for w in words:
    d = sym_sub(sym_expanded(sym_G(w, c0)), sym_expanded(sym6(sp.Integer(0), w, c0)))
    ok, rk = decomposable(d, [sp.Integer(2), sp.Integer(3)])
    print(f"word {w}: mine-cor6 decomposable: {ok}   (diff terms: {len(d)})")
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
