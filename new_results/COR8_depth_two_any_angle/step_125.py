
import sys, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
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
nonconst = [(s, rv, omv) for s, rv, omv in sa if rv != Z4]
pairs = []
for i, (s1, rv1, omv1) in enumerate(nonconst):
    for j, (s2, rv2, omv2) in enumerate(nonconst):
        v = factor_vec(sp.expand(1 - sp.sympify(s1) * sp.sympify(s2)))
        if v is not None:
            pairs.append((s1, s2, rv1, rv2, omv1, omv2, v))
print(len(pairs), "pairs")
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
# target: the [2]-residual of delta22(Omega-DT)
target = {((0,2),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,2),(1,3)): {0: Fr(1), -1: Fr(1)},
          ((0,3),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,3),(1,3)): {0: Fr(1), -1: Fr(1)}}
rows = []
for (s1, s2, rv1, rv2, omv1, omv2, rom) in pairs:
    for kind, fn in ((22, li22_sym), (31, li31_sym)):
        S = dropz(fn(rv1, rv2, omv1, omv2, rom))
        d = delta22(S, 4)
        d2 = {k: v for k, v in d.items() if k in target or any(i == 0 for kk in k for i in kk)}
        if d:
            rows.append((f"Li{kind}({s1},{s2})", d))
print(len(rows), "nonzero delta22")
# those with support exactly in target's 4 components
exact = [(nm, d) for nm, d in rows if set(d.keys()) == set(target.keys())]
print(len(exact), "with exact support")
for nm, d in exact[:20]:
    print(nm, {k: {e: str(v) for e,v in sorted(p.items())} for k,p in sorted(d.items())})
pickle.dump([(nm, {str(k): {str(e): str(v) for e,v in p.items()} for k,p in d.items()}) for nm, d in rows], open('pairs_delta22_2kept.pkl','wb'))
