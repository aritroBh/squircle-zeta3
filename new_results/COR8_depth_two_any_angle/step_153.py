
import sys, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
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
Rm = ladd(pv, lneg(L2))
for nm, S in [('Li22(1-c,1/2)', dropz(li22_sym(mv, lneg(L2), cv, lneg(L2), Rm))),
              ('Li31(1-c,1/2)', dropz(li31_sym(mv, lneg(L2), cv, lneg(L2), Rm)))]:
    v = to256(S)
    print(nm)
    for w in ['m2cm','mcm2','mcmp','mcpm','2mcm','2pcm','pcm2','pcmp','pcpm','pmcm','p2cm','2cm2','c22m','c2m2','cm22']:
        val = v[widx(w)]
        if val:
            print("   ", w, val)
