t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
from functools import lru_cache
import sympy as sp, json, re
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

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = sorted(set(parse_term(s)[1] for s, c in data if len(parse_term(s)[1]) == 4))
x = sp.Symbol('x')
nbad = 0
for w in W4:
    s1 = sym_G(w, x)
    s2 = sym6(sp.Integer(0), w, x)
    if sym_sub(s1, s2) != {}:
        nbad += 1
        print("DIFF at word", w)
        print("   mine :", s1)
        print("   cor6 :", s2)
print("generic x: words differing:", nbad, "of", len(W4))
nbad = 0
c0 = sp.Rational(1, 2)
for w in W4:
    s1 = sym_G(w, c0)
    s2 = sym6(sp.Integer(0), w, c0)
    if sym_sub(s1, s2) != {}:
        nbad += 1
        print("DIFF at 1/2 word", w)
        print("   mine :", s1)
        print("   cor6 :", s2)
print("at 1/2: words differing:", nbad, "of", len(W4))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
