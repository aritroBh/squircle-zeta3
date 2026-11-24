
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
Rm = (Fr(-1), Fr(0), Fr(0), Fr(1)); Rp = (Fr(-1), Fr(0), Fr(1), Fr(0))
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
func = {'222m': 1, '22m2': -1, '22mm': -1, '22mp': -1, 'm22m': 1, 'p22m': 1}
def pairing(S):
    v = to256(S)
    return sum(c * v[widx(w)] for w, c in func.items())
S22m = dropz(li22_sym(lneg(L2), mv, lneg(L2), cv, Rm))
S22p = dropz(li22_sym(lneg(L2), pv, lneg(L2), cv, Rp))
print("pairing Li22(1/2,1-c):", pairing(S22m))
print("pairing Li22(1/2,1+c):", pairing(S22p))
# check word indices
print("word '222m' =", widx('222m'), "=", bwords[widx('222m')])
print("word 'm22m' =", widx('m22m'), "=", bwords[widx('m22m')])
# print actual content of S22m at those words
v = to256(S22m)
for w, c in func.items():
    print(f"  {w}: {v[widx(w)]}")
