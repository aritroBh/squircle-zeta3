"""step_02.py -- delta22 matching: enumerate (g,h) S-unit pairs, compute delta22 of
Li31(g,h) and Li22(g,h), solve sum b_j(c) delta_j = delta22(Omega), minimal terms."""
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *

c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"

# ---------- S-unit factorisation (as in step13_library) ----------
def factor_vec(rat):
    num, den = sp.fraction(sp.together(rat))
    vec = [0, 0, 0, 0]  # [2, c, 1-c, 1+c]
    for poly, sgn in [(sp.expand(num), +1), (sp.expand(den), -1)]:
        if poly == 0:
            return None
        fac = sp.factor_list(poly)
        base = sp.Rational(fac[0])
        b = abs(base)
        if b != 1:
            bn, bd = b.numerator, b.denominator
            if bn & (bn - 1) or bd & (bd - 1):
                return None
            vec[0] += sgn * (bn.bit_length() - bd.bit_length())
        for f, e in fac[1]:
            f = sp.expand(f)
            if f == c:
                vec[1] += sgn * e
            elif f == c - 1:
                vec[2] -= sgn * e
            elif f == c + 1:
                vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f))
                fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1):
                    return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else:
                return None
    return tuple(Fr(v) for v in vec)

sa = json.load(open(C4 + "sunit_args.json"))
args = []  # (expr, vec_r3, vec_om3) with the [2]-component dropped
for v_r, v_om, s in sa:
    ex = sp.sympify(s)
    args.append((ex, tuple(Fr(x) for x in v_r[1:]), tuple(Fr(x) for x in v_om[1:])))
print(f"{len(args)} S-unit args")

Z3 = (Fr(0),) * 3
# pairs with 1 - g*h factoring
pairs = []
for i, (g, gv, gom) in enumerate(args):
    for j, (h, hv, hom) in enumerate(args):
        v = factor_vec(sp.expand(1 - g * h))
        if v is not None:
            pairs.append((i, j, gv, hv, gom, hom, v[1:]))
print(f"{len(pairs)} pairs with 1-gh an S-unit")

def delta_of(gv, hv, gom, hom, romgh, kind):
    X, Y, P, Q, R = gv, hv, gom, hom, romgh
    S = li31_sym(X, Y, P, Q, R) if kind == 31 else li22_sym(X, Y, P, Q, R)
    return delta22(S, 3)

target = {((0, 1), (1, 2)): {0: Fr(1), -1: Fr(-1)},
          ((0, 2), (1, 2)): {0: Fr(-1), -1: Fr(-1)}}

# compute all deltas, collect distinct ones with support inside the 3 components
rows = []
for (i, j, gv, hv, gom, hom, r) in pairs:
    for kind in (31, 22):
        d = delta_of(gv, hv, gom, hom, r, kind)
        if d:
            rows.append((i, j, kind, d))
print(f"{len(rows)} nonzero delta22 values")
from collections import Counter
supp = Counter(tuple(sorted(d.keys())) for _, _, _, d in rows)
for s, n in supp.most_common():
    print(f"support {s}: {n}")
json.dump(rows, open("pairs_delta22.json", "w"))
