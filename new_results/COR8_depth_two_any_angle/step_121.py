
import sys, json, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
Z4 = (Fr(0),)*4   # basis ([2],[c],[1-c],[1+c])
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)

def factor_vec(rat):
    num, den = sp.fraction(sp.together(rat))
    vec = [0, 0, 0, 0]
    for poly, sgn in [(sp.expand(num), +1), (sp.expand(den), -1)]:
        if poly == 0: return None
        fac = sp.factor_list(poly)
        base = sp.Rational(fac[0]); b = abs(base)
        if b != 1:
            bn, bd = b.numerator, b.denominator
            if bn & (bn - 1) or bd & (bd - 1): return None
            vec[0] += sgn * (bn.bit_length() - bd.bit_length())
        for f, e in fac[1]:
            f = sp.expand(f)
            if f == c: vec[1] += sgn * e
            elif f == c - 1: vec[2] += sgn * e
            elif f == c + 1: vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f)); fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1): return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else: return None
    return tuple(Fr(v) for v in vec)

args = {}
for a in range(-2, 3):
    for b in range(-3, 4):
        for d in range(-3, 4):
            for e in range(-3, 4):
                for sgn in [1, -1]:
                    r = sgn * sp.Integer(2)**a * c**b * (1 - c)**d * (1 + c)**e
                    v_r = (Fr(a), Fr(b), Fr(d), Fr(e))
                    v_om = factor_vec(sp.expand(1 - r))
                    if v_om is not None:
                        args[(v_r, v_om)] = r
sa = [(sp.sstr(r), v_r, v_om) for (v_r, v_om), r in sorted(args.items(), key=lambda kv: str(kv[1]))]
polefree = [(s, rv, omv) for s, rv, omv in sa if rv != Z4 and sp.sympify(s).subs(c, 0).is_finite]
print(len(polefree), "pole-free nonconstant args")
# F symbol with [2] kept
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
# DT symbol with [2]: S(Li22(1/2,y)) = R Y Y Y - R Y Q Y, R=[1-y/2], Q=[1-y]
# y=1-c: R=[(1+c)/2]=pv-L2, Y=mv, Q=cv ; y=1+c: R=[(1-c)/2]=mv-L2, Y=pv, Q=cv
def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
Rm = (Fr(-1), Fr(0), Fr(0), Fr(1))
Rp = (Fr(-1), Fr(0), Fr(1), Fr(0))
S22m = drop_zero(li22_sym(L2, mv, L2, cv, Rm))
S22p = drop_zero(li22_sym(L2, pv, L2, cv, Rp))
print("S(Li22(1/2,1-c)) with [2]:", len(S22m), "terms;", len(S22p), "terms")
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
print("S(Omega-DT) with [2]:", len(R), "words")
pickle.dump({'OM': OM, 'DT': DT, 'R': R, 'polefree': polefree}, open('symbols_2kept.pkl','wb'))
