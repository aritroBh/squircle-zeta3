
import sys, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
Rm = ladd(pv, lneg(L2))
Rp = ladd(mv, lneg(L2))
for nm, S in [('Li31(1-c,1/2)', dropz(li31_sym(mv, lneg(L2), cv, lneg(L2), Rm))),
              ('Li31(1+c,1/2)', dropz(li31_sym(pv, lneg(L2), mv, lneg(L2), Rp))),
              ('Li31(1-c,2)', dropz(li31_sym(mv, L2, cv, lneg(L2), lneg(cv)))),
              ('Li31(1+c,2)', dropz(li31_sym(pv, L2, mv, lneg(L2), lneg(cv))))]:
    print(nm)
    for w, cf in sorted(S.items(), key=str):
        def ws(w):
            out = []
            for l in w:
                if sum(map(abs, l)) == 1:
                    out.append('2cmp'[next(j for j, x in enumerate(l) if x)])
                else:
                    out.append(str(tuple(int(x) for x in l)))
            return ''.join(out)
        print("   ", ws(w), {e: str(v) for e, v in cf.items()})
