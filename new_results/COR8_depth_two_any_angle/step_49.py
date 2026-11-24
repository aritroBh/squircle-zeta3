
import sys, json, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
import itertools
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
Z3 = (Fr(0),)*3
def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}

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
            elif f == c - 1: vec[2] -= sgn * e
            elif f == c + 1: vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f)); fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1): return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else: return None
    return tuple(Fr(v) for v in vec)

sa = json.load(open(C4 + "sunit_args.json"))
args = []
for v_r, v_om, s in sa:
    args.append((sp.sympify(s), tuple(Fr(x) for x in v_r[1:]), tuple(Fr(x) for x in v_om[1:])))
pairs = []
for i, (g, gv, gom) in enumerate(args):
    for j, (h, hv, hom) in enumerate(args):
        v = factor_vec(sp.expand(1 - g * h))
        if v is not None:
            pairs.append((sp.sstr(g), sp.sstr(h), gv, hv, gom, hom, v[1:]))
print(len(pairs), "pairs")
bwords = list(itertools.product(range(3), repeat=4))
bidx = {w:i for i,w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(3) if l[j]]
def to81(S):
    v = [Fr(0)]*81
    for w, cf in S.items():
        coeff = cf.get(0, Fr(0))
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for pp, c0 in terms:
            v[bidx[pp]] += coeff*c0
    return v
dt_syms = {}
for (gs, hs, gv, hv, gom, hom, r) in pairs:
    for kind, fn in ((22, li22_sym), (31, li31_sym)):
        key = f"Li{kind}({gs},{hs})"
        S = drop_zero(fn(gv, hv, gom, hom, r))
        v = to81(S)
        if any(v):
            dt_syms[key] = v
print(len(dt_syms), "depth-two candidate symbols")
# dedupe identical symbols
seen = {}
dt_list = []
for k, v in dt_syms.items():
    t = tuple(v)
    if t not in seen:
        seen[t] = k; dt_list.append((k, v))
print(len(dt_list), "distinct")
pickle.dump(dt_list, open('dt_syms.pkl','wb'))
open('progress.log','a').write(f'step_03b: {len(dt_list)} distinct depth-two candidate symbols\n')
