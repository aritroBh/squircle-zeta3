
import sys, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
c = sp.Symbol('c')
Z4 = (Fr(0),)*4
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
for a in range(-3, 4):
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
nonconst = [(s, rv, omv) for s, rv, omv in sa if rv != Z4]
print(len(sa), "args", len(nonconst), "nonconst")
pairs = []
for i, (s1, rv1, omv1) in enumerate(nonconst):
    for j, (s2, rv2, omv2) in enumerate(nonconst):
        v = factor_vec(sp.expand(1 - sp.sympify(s1) * sp.sympify(s2)))
        if v is not None:
            pairs.append((s1, s2, rv1, rv2, omv1, omv2, v))
print(len(pairs), "pairs")
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
hits = []
for (s1, s2, rv1, rv2, omv1, omv2, rom) in pairs:
    for kind, fn in ((22, li22_sym), (31, li31_sym)):
        S = dropz(fn(rv1, rv2, omv1, omv2, rom))
        v = to256(S)
        if not any(v): continue
        s1_ = sum(co * v[widx(w)] for w, co in F1.items())
        s2_ = sum(co * v[widx(w)] for w, co in F2.items())
        if s1_ or s2_:
            hits.append((f"Li{kind}({s1},{s2})", s1_, s2_))
for nm, a, b in hits:
    print(nm, "F1:", a, " F2:", b)
