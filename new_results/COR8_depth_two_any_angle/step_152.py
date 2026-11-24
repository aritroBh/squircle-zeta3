
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
# Li22(1-c, 1/2): X=mv, P=cv, Y=-[2], Q=-[2], R=pv-[2]
Rm = ladd(pv, lneg(L2))
S = dropz(li22_sym(mv, lneg(L2), cv, lneg(L2), Rm))
v = to256(S)
F2 = {'2mcm': 1, '2pcm': 1, 'm2cm': -1, 'pcm2': -1, 'pcmp': -1, 'pcpm': 1}
F1 = {'m2cm': 1, 'mcm2': -1, 'mcmp': -1, 'mcpm': 1}
print("Li22(1-c,1/2) F2 pairing:", sum(co*v[widx(w)] for w, co in F2.items()))
print("Li22(1-c,1/2) F1 pairing:", sum(co*v[widx(w)] for w, co in F1.items()))
# and Li31(1-c, 1/2)
S31 = dropz(li31_sym(mv, lneg(L2), cv, lneg(L2), Rm))
v31 = to256(S31)
print("Li31(1-c,1/2) F2:", sum(co*v31[widx(w)] for w, co in F2.items()), " F1:", sum(co*v31[widx(w)] for w, co in F1.items()))
# Li22(1+c, 1/2): X=pv, P=mv, R=mv-[2]
Rp = ladd(mv, lneg(L2))
Sp = dropz(li22_sym(pv, lneg(L2), mv, lneg(L2), Rp))
vp = to256(Sp)
print("Li22(1+c,1/2) F2:", sum(co*vp[widx(w)] for w, co in F2.items()), " F1:", sum(co*vp[widx(w)] for w, co in F1.items()))
S31p = dropz(li31_sym(pv, lneg(L2), mv, lneg(L2), Rp))
v31p = to256(S31p)
print("Li31(1+c,1/2) F2:", sum(co*v31p[widx(w)] for w, co in F2.items()), " F1:", sum(co*v31p[widx(w)] for w, co in F1.items()))
