
import sys, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
def lsub(u, v): return ladd(u, lneg(v))
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
bwords = list(itertools.product(range(4), repeat=4))
bidx = {w: i for i, w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(4) if l[j]]
def to256(S, e=0):
    v = [Fr(0)]*256
    for w, cf in S.items():
        coeff = cf.get(e)
        if coeff is None: continue
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for pp, c0 in terms:
            v[bidx[pp]] += coeff*c0
    return v
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
F1 = {'m2cm': 1, 'mcm2': -1, 'mcmp': -1, 'mcpm': 1}
F2 = {'2mcm': 1, '2pcm': 1, 'm2cm': -1, 'pcm2': -1, 'pcmp': -1, 'pcpm': 1}
def pairing(S, F):
    v = to256(S)
    return sum(co*v[widx(w)] for w, co in F.items())
cands = [
    ('((1-c)/2, 1/(1+c))', (lsub(mv, L2), lneg(pv), lsub(pv, L2), lsub(cv, pv), lsub(cv, pv))),
    ('((1-c)/2, 2/(1+c))', (lsub(mv, L2), lsub(L2, pv), lsub(pv, L2), lsub(mv, pv), ladd(L2, lsub(cv, pv)))),
    ('((1+c)/2, 2c/(1+c))', (lsub(pv, L2), ladd(L2, lsub(cv, pv)), lsub(mv, L2), lsub(mv, pv), cv)),
    ('((1-c)/2, 2c/(1-c))', (lsub(mv, L2), ladd(L2, lsub(cv, mv)), lsub(pv, L2), lsub(pv, mv), cv)),
]
for nm, (X, Y, P, Q, R) in cands:
    for kind, fn in ((22, li22_sym), (31, li31_sym)):
        S = dropz(fn(X, Y, P, Q, R))
        print(f"{nm} Li{kind}: F1={pairing(S, F1)}  F2={pairing(S, F2)}")
