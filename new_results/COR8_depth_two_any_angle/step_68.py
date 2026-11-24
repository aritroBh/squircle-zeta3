
import json
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
def factor_vec(rat):
    """CORRECT: [c-1] = [1-c] (sign ignored)."""
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
            elif f == c - 1: vec[2] += sgn * e   # FIXED: [c-1]=[1-c]
            elif f == c + 1: vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f)); fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1): return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else: return None
    return tuple(Fr(v) for v in vec)

sa = json.load(open("<project>/research/08_theory/COR4_M3alpha_classical/sunit_args.json"))
bad = 0
for v_r, v_om, s in sa:
    ex = sp.sympify(s)
    rv = factor_vec(ex); omv = factor_vec(1 - ex)
    if tuple(Fr(x) for x in v_r) != rv or tuple(Fr(x) for x in v_om) != omv:
        print("WRONG in sunit_args:", s, "stored:", v_r, v_om, "correct:", rv, omv)
        bad += 1
print("wrong entries in sunit_args.json:", bad, "of", len(sa))
