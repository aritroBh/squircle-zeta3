
import sys, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
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
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
# the violated functional
func_words = ['222m','22m2','22mm','22mp','m22m','p22m']
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
# candidate: Li31(1/2, 1-c) and Li31(1/2, 1+c), Li22 same
Rm = (Fr(-1), Fr(0), Fr(0), Fr(1)); Rp = (Fr(-1), Fr(0), Fr(1), Fr(0))
for nm, S in [
    ('Li31(1/2,1-c)', dropz(li31_sym(lneg(L2), mv, lneg(L2), cv, Rm))),
    ('Li31(1/2,1+c)', dropz(li31_sym(lneg(L2), pv, lneg(L2), cv, Rp))),
    ('Li22(1/2,1-c)', dropz(li22_sym(lneg(L2), mv, lneg(L2), cv, Rm))),
    ('Li22(1/2,1+c)', dropz(li22_sym(lneg(L2), pv, lneg(L2), cv, Rp))),
]:
    v = to256(S)
    pairing = [v[widx(w)] for w in func_words]
    print(nm, "pairing with violated functional:", [str(x) for x in pairing])
# also Li22 with y=2c/(1-c), y=2c/(1+c) (DGR): letters
# x=1/2: X=-[2], P=-[2]; y=2c/(1-c): Y=[2]+[c]-[m], Q=[1-y]=[(1+c)/(1-c)]=pv-mv, xy=c/(c-1): R=[1-xy]=[( -1)/(c-1)]= -[m]... = lneg(mv)
Yd = ladd(L2, ladd(cv, lneg(mv))); Qd = lsub(pv, mv); Rd = lneg(mv)
Yp = ladd(L2, ladd(cv, lneg(pv))); Qp = lsub(mv, pv); Rp2 = lneg(pv)
for nm, S in [
    ('Li22(1/2,2c/(c-1))', dropz(li22_sym(lneg(L2), Yd, lneg(L2), Qd, Rd))),
    ('Li22(1/2,2c/(c+1))', dropz(li22_sym(lneg(L2), Yp, lneg(L2), Qp, Rp2))),
]:
    v = to256(S)
    pairing = [v[widx(w)] for w in func_words]
    print(nm, "pairing:", [str(x) for x in pairing])
