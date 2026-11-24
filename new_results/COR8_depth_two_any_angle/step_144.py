
import sys, pickle, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
# residual to match:
target = {((0,2),(1,2)): {0: Fr(-1), -1: Fr(-1)},
          ((0,2),(1,3)): {0: Fr(-1), -1: Fr(-1)},
          ((0,3),(1,2)): {0: Fr(-1), -1: Fr(1)},
          ((0,3),(1,3)): {0: Fr(1), -1: Fr(1)}}
# build delta22 of log2-terms for all args
import json
sa = json.load(open("<project>/research/08_theory/COR4_M3alpha_classical/sunit_args.json"))
c = sp.Symbol('c')
args = [s for _,_,s in sa]
logs = {'c': cv, '1-c': mv, '1+c': pv}
rows = []
import sympy as sp
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
for s in args:
    ex = sp.sympify(s)
    rv = factor_vec(ex); omv = factor_vec(1 - ex)
    if rv is None or omv is None: continue
    S3 = li_sym(3, rv, omv); S2 = li_sym(2, rv, omv)
    for nm, S in [(f"l2*Li3({s})", sshuffle(log_sym(L2), S3)),
                  (f"l2^2*Li2({s})", sshuffle(sshuffle(log_sym(L2), log_sym(L2)), S2))]:
        rows.append((nm, delta22(dropz(S), 4)))
    for hn, hv in logs.items():
        rows.append((f"l2*Li2({s})*log({hn})", delta22(dropz(sshuffle(log_sym(L2), sshuffle(S2, log_sym(hv)))), 4)))
import itertools
for combo in itertools.combinations_with_replacement(list(logs), 3):
    S = log_sym(L2)
    for h in combo:
        S = sshuffle(S, log_sym(logs[h]))
    rows.append(("l2*" + "*".join(f"log({h})" for h in combo), delta22(dropz(S), 4)))
for combo in itertools.combinations_with_replacement(list(logs), 2):
    S = sshuffle(log_sym(L2), log_sym(L2))
    for h in combo:
        S = sshuffle(S, log_sym(logs[h]))
    rows.append(("l2^2*" + "*".join(f"log({h})" for h in combo), delta22(dropz(S), 4)))
for hn, hv in logs.items():
    S = sshuffle(sshuffle(log_sym(L2), log_sym(L2)), log_sym(L2))
    rows.append((f"l2^3*log({hn})", delta22(dropz(sshuffle(S, log_sym(hv))), 4)))
rows = [(nm, d) for nm, d in rows if d]
print(len(rows), "log2-terms with nonzero delta22")
for nm, d in rows[:15]:
    print(nm, {k: {e: str(v) for e, v in sorted(p.items())} for k, p in sorted(d.items())})
pickle.dump([(nm, {str(k): {str(e): str(v) for e, v in p.items()} for k, p in d.items()}) for nm, d in rows], open('log2_delta22.pkl','wb'))
