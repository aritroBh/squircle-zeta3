t = r'''
# step_12.py -- compare my recursion vs COR6's documented recursion
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
from functools import lru_cache
import sympy as sp

# COR6 formula (from their cor6_sym.py docstring):
#   S(J(a1..an;x)) = sum_j [l_j] (x) S(J(omit a_j)),  l_j = (a_{j+1}-a_j)/(a_{j-1}-a_j),
#   a_0 := 0, a_{n+1} := x. Zero numerator/denominator factor dropped; both vanish -> term 0.
# They prepend l_j: tensor order has l_j FIRST.
# eps-regularization: perturb a_i -> a_i + eps^i (same as mine).
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
        if P[0] == 'c' and Q[0] == 'c':
            raise AssertionError("endpoints collide")
        return EPS
    return canon6(d)

@lru_cache(maxsize=None)
def sym6(a0, letters, aN):
    n = len(letters)
    A0 = ('c', a0); AN = ('c', aN)
    LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
    if n == 0:
        return {(): 1}
    out = {}
    for j in range(n):
        p = A0 if j == 0 else LL[j - 1]
        q = AN if j == n - 1 else LL[j + 1]
        num = dlog6(q, LL[j])   # a_{j+1} - a_j
        den = dlog6(p, LL[j])   # a_{j-1} - a_j
        subs = letters[:j] + letters[j+1:]
        Ssub = sym6(a0, subs, aN)
        lv = []
        if not (num is EPS or iszero6(num)): lv.append(+1)
        if not (den is EPS or iszero6(den)): lv.append(-1)
        for sgn, ent in ((+1, num), (-1, den)):
            if ent is EPS or iszero6(ent):
                continue
            for w2, c2 in Ssub.items():
                out[(ent,) + w2] = out.get((ent,) + w2, 0) + sgn * c2
    return {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}

x = sp.Symbol('x')
print("their S(G(0,1;x)) =", sym6(sp.Integer(0), (sp.Integer(0), sp.Integer(1)), x))
print("mine  S(G(0,1;x)) =", sym_G((0,1), x))
print("their S(G(0,0,0,1;x)) =", sym6(sp.Integer(0), (0,0,0,1), x))
print("mine  S(G(0,0,0,1;x)) =", sym_G((0,0,0,1), x))
print("their S(Li_2(x)) (=-G(0,1)) =", sym6(sp.Integer(0), (sp.Integer(0), sp.Integer(1)), x))
print("their S(G(1,0;x)) =", sym6(sp.Integer(0), (sp.Integer(1), sp.Integer(0)), x))
print("mine  S(G(1,0;x)) =", sym_G((1,0), x))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
