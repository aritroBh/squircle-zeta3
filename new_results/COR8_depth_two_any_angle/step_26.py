
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"

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
            pairs.append((i, j, gv, hv, gom, hom, v[1:]))

W1 = ((0,1),(1,2)); W2 = ((0,2),(1,2)); W12 = ((0,1),(0,2))
def delta_of(gv, hv, gom, hom, r, kind):
    S = li31_sym(gv, hv, gom, hom, r) if kind == 31 else li22_sym(gv, hv, gom, hom, r)
    return delta22(S, 3)

rows = []
for (i, j, gv, hv, gom, hom, r) in pairs:
    for kind in (31, 22):
        d = delta_of(gv, hv, gom, hom, r, kind)
        if d and set(d.keys()) <= {W1, W2}:
            rows.append((i, j, kind, d))
print(len(rows), "pairs with support in {W1,W2}")
# print distinct values
from collections import defaultdict
byval = defaultdict(list)
for (i,j,k,d) in rows:
    key = tuple(sorted((kk, tuple(sorted(vv.items()))) for kk, vv in d.items()))
    byval[key].append((i,j,k))
for key, occ in sorted(byval.items(), key=lambda kv: str(kv[0])):
    d = {k: dict(v) for k, v in key}
    s = {('W1' if k==W1 else 'W2'): {e: str(v) for e,v in sorted(vv.items())} for k,vv in d.items()}
    ex = occ[:3]
    exs = [f"Li{k}({sp.sstr(args[i][0])},{sp.sstr(args[j][0])})" for i,j,k in ex]
    print(s, f"  x{len(occ)}", exs)
