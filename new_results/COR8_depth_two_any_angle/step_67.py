
import pickle, sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
import itertools
c = sp.Symbol('c')
Z3 = (Fr(0),)*3
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
    return tuple(Fr(v) for v in vec[1:])
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv = D['li4v'], D['prodv']
li4n = list(li4v.keys()); prodn = list(prodv.keys())
# check factor_vec on all Li4 arg names vs the stored sunit vectors
sa = json.load(open("<project>/research/08_theory/COR4_M3alpha_classical/sunit_args.json"))
sunit = {s: (tuple(Fr(x) for x in v_r[1:]), tuple(Fr(x) for x in v_om[1:])) for v_r, v_om, s in sa}
bad = 0
for nm in li4n:
    rv = factor_vec(sp.sympify(nm)); omv = factor_vec(1 - sp.sympify(nm))
    if (rv, omv) != sunit[nm]:
        print("MISMATCH", nm, rv, omv, "vs", sunit[nm]); bad += 1
print("bad Li4 args:", bad)
