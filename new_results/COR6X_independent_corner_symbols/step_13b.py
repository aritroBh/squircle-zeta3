
# step_13b.py -- number identity with the coproduct recursion C
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
    n = len(letters)
    A0 = ('c', a0); AN = ('c', aN)
    LL = tuple(('m', a, i + 1) for i, a in enumerate(letters))
    if n == 0: return {(): 1}
    out = {}
    for i in range(n):
        ai = LL[i]
        e1 = dl(AN, ai); e0 = dl(A0, ai)
        SL = symC(a0, tuple(l[1] for l in LL[:i]), ai[1]) if i > 0 else {(): 1}
        SR = symC(ai[1], tuple(l[1] for l in LL[i+1:]), aN) if i < n - 1 else {(): 1}
        Sseg = sym_shuffle(SL, SR)
        for sgn, ent in ((+1, e1), (-1, e0)):
            if ent is EPS or isz(ent): continue
            for w2, c2 in Sseg.items():
                out[(ent,) + w2] = out.get((ent,) + w2, 0) + sgn * c2
    return {k: v for k, v in out.items() if v and not any(t is EPS for t in k)}
def sC_G(letters, pt): return symC(sp.Integer(0), letters, pt)
def sC_H(word, pt):
    sgn = (-1) ** sum(1 for a in word if a == 1)
    return {k: sgn * v for k, v in sC_G(word, pt).items()}

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
c0 = sp.Rational(1, 2)
S_om = {}
for w, cf in TERMS:
    S_om = sym_add(S_om, sym_scale(sC_H(w, c0), cf))
S_om = sym_expanded(S_om)
def sC_Li4(z): return sym_scale(sC_G((0,0,0,1), z), -1)
S_five = sym_expanded(sym_add(
    sym_scale(sC_Li4(-2), sp.Rational(-21)),
    sym_scale(sC_Li4(sp.Rational(1,3)), sp.Rational(-65,8)),
    sym_scale(sC_Li4(sp.Rational(2,3)), sp.Rational(13,2)),
    sym_scale(sC_Li4(sp.Rational(1,2)), sp.Rational(137,6)),
    sym_scale(sC_Li4(sp.Rational(1,9)), sp.Rational(39,64))))
S_log2 = {(sp.Integer(2),): 1}; S_log3 = {(sp.Integer(3),): 1}
S_Li3_13 = sym_expanded(sym_scale(sC_G((0,0,1), sp.Rational(1,3)), -1))
S_Li3_m12 = sym_expanded(sym_scale(sC_G((0,0,1), sp.Rational(-1,2)), -1))
S_Li2_13 = sym_expanded(sym_scale(sC_G((0,1), sp.Rational(1,3)), -1))
def sh(*args):
    r = args[0]
    for a in args[1:]: r = sym_shuffle(r, a)
    return r
S_P = sym_add(
    sym_scale(sh(S_Li3_13, S_log2), sp.Rational(15,2)),
    sym_scale(sh(S_Li3_m12, S_log2), sp.Rational(37,2)),
    sym_scale(sh(S_Li3_m12, S_log3), sp.Rational(-5,4)),
    sh(S_Li2_13, S_Li2_13),
    sym_scale(sh(S_Li2_13, S_log2, S_log2), sp.Rational(-23,4)),
    sym_scale(sh(S_Li2_13, S_log2, S_log3), sp.Rational(7,4)),
    sh(S_Li2_13, S_log3, S_log3),
    sym_scale(sh(S_log2, S_log2, S_log2, S_log2), sp.Rational(-211,144)),
    sym_scale(sh(S_log2, S_log2, S_log2, S_log3), sp.Rational(101,24)),
    sym_scale(sh(S_log2, S_log2, S_log3, S_log3), sp.Rational(-15,8)),
    sym_scale(sh(S_log2, S_log3, S_log3, S_log3), sp.Rational(-35,24)),
    sym_scale(sh(S_log3, S_log3, S_log3, S_log3), sp.Rational(113,192)))
S_P = sym_expanded(S_P)
T = sym_sub(S_om, sym_add(S_five, S_P))
print("T_C = S_C(Omega(1/2)) - S_C(five) - S_C(P): terms =", len(T))
for k, v in sorted(T.items(), key=str): print("   ", k, v)
