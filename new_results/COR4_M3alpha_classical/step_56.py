cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
step = '''"""step13_library.py -- enumerate S-unit arguments r(c): r and 1-r both factor
over {2, c, 1-c, 1+c}; report the list with letter vectors."""
import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
import sympy as sp

c = sp.Symbol('c')
# letter basis: 0:[2], 1:[c], 2:[1-c], 3:[1+c]
def factor_vec(rat):
    """rat: sympy rational function of c. If rat = +-2^a c^b (1-c)^d (1+c)^e,
    return (a,b,d,e) else None."""
    num, den = sp.fraction(sp.expand(rat))
    if not num.is_polynomial(c) or not den.is_polynomial(c):
        return None
    vec = [0, 0, 0, 0]
    for poly, sgn in [(num, +1), (den, -1)]:
        poly = sp.Poly(poly, c)
        coefs = poly.all_coeffs()
        # constant factor
        cont = sp.Rational(coefs[0])
        # factor the polynomial
        fac = sp.factor_list(sp.expand(poly.as_expr()))
        base, pairs = fac
        base = sp.Rational(base)
        # base must be +- power of 2
        b = abs(sp.Rational(base))
        if b != 1:
            if b.denominator != 1 and b.numerator != 1:
                pass
            bn, bd = b.numerator, b.denominator
            if bn & (bn - 1) or bd & (bd - 1):
                return None
            import math
            vec[0] += sgn * (bn.bit_length() - bd.bit_length())
        for f, e in pairs:
            f = sp.expand(f)
            if f == c:
                vec[1] += sgn * e
            elif f == c - 1:
                vec[2] += sgn * e * -1   # (c-1)^e = (-1)^e (1-c)^e
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

args = {}
for a in range(-2, 3):
    for b in range(-3, 4):
        for d in range(-3, 4):
            for e in range(-3, 4):
                for sgn in [1, -1]:
                    r = sgn * 2**a * c**b * (1 - c)**d * (1 + c)**e
                    v_r = (Fr(a), Fr(b), Fr(d), Fr(e))
                    om = sp.expand(1 - r)
                    v_om = factor_vec(om)
                    if v_om is not None:
                        args[(v_r, v_om)] = r
print(f"{len(args)} S-unit arguments (r, 1-r) over {{2,c,1-c,1+c}}")
for (v_r, v_om), r in sorted(args.items(), key=lambda kv: str(kv[1])):
    print(f"  r = {str(r):>28}   [r]={v_r}  [1-r]={v_om}")
import json
json.dump([[list(map(str, v_r)), list(map(str, v_om)), str(r)] for (v_r, v_om), r in args.items()],
          open("<project>/research/08_theory/COR4_M3alpha_classical/sunit_args.json", "w"), indent=1)
'''
open(cor+"step13_library.py","w").write(step)
import sys, subprocess
r = subprocess.run([sys.executable, cor+"step13_library.py"], capture_output=True, text=True, timeout=300)
print(r.stdout)
print(r.stderr[-800:] if r.returncode else "exit 0")
