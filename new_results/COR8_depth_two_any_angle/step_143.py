
import sys, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
import json
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
FS = json.load(open(C4 + "F_symbol.json"))
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            digs.append(int(s[i:i+2])); i += 2
        else:
            digs.append(int(s[i])); i += 1
    assert digs[0] == 0
    return (Fr(digs[0]), Fr(digs[1]), Fr(digs[2]), Fr(digs[3]))
OM = {}
for key, val in FS.items():
    e = int(key.split(" ")[0][2:])
    word = tuple(parse_letter(b) for b in key.split(" ")[1].split("|"))
    OM[word] = padd(OM.get(word, {}), {e: Fr(val)})
def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
Rm = (Fr(-1), Fr(0), Fr(0), Fr(1))
Rp = (Fr(-1), Fr(0), Fr(1), Fr(0))
h2 = lneg(L2)
S22m = drop_zero(li22_sym(h2, mv, h2, cv, Rm))
S22p = drop_zero(li22_sym(h2, pv, h2, cv, Rp))
DT = {}
for w, cf in S22m.items():
    DT[w] = pmul(cf, {-1: Fr(1), 0: Fr(-1)})
for w, cf in S22p.items():
    DT[w] = padd(DT.get(w, {}), pmul(cf, {-1: Fr(-1), 0: Fr(-1)}))
DT = {w: cf for w, cf in DT.items() if cf}
R = dict(OM)
for w, cf in DT.items():
    R[w] = padd(R.get(w, {}), pmul(cf, {0: Fr(-1)}))
R = {w: cf for w, cf in R.items() if cf}
print("S(Omega-DT) [2]-kept:", len(R), "words")
dR = delta22(R, 4)
print("delta22(Omega-DT):", {k: {e: str(v) for e, v in sorted(p.items())} for k, p in sorted(dR.items())})
pickle.dump({'OM': OM, 'DT': DT, 'R': R, 'polefree': D['polefree'] if 'polefree' in D else None}, open('symbols_2kept.pkl','wb'))
# check pairing with the functional
import itertools
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
Rv = {E: to256(R, E) for E in (-1, 0, 1)}
func = {'222m': 1, '22m2': -1, '22mm': -1, '22mp': -1, 'm22m': 1, 'p22m': 1}
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
for E in (-1, 0, 1):
    print(f"f . R_{E} =", sum(c * Rv[E][widx(w)] for w, c in func.items()))
pickle.dump(Rv, open('Rv_2kept.pkl','wb'))
