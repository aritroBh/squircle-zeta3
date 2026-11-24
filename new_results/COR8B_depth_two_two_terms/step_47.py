
import sys, pickle, itertools
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr
import sympy as sp
c = sp.symbols('c')

# rebuild the 38-arg S-unit list (unsigned values + signs) quickly from pickle of pairs? 
# Simpler: reuse pairs_delta_all_v2 args by parsing strings is fragile; regenerate with smaller grid (deg<=2 by Mason)
def factor_vec(expr):
    expr = sp.simplify(expr)
    num, den = sp.fraction(sp.cancel(expr))
    vec = [0,0,0,0]
    for poly, sgn in [(num,1),(den,-1)]:
        poly = sp.Poly(poly, c)
        content, factors = sp.factor_list(poly)
        cont = sp.Rational(content);  cont = -cont if cont<0 else cont
        for prime, e in sp.factorint(int(cont.p)).items():
            if prime != 2: return None
            vec[0] += sgn*e
        for prime, e in sp.factorint(int(cont.q)).items():
            if prime != 2: return None
            vec[0] -= sgn*e
        for f, e in factors:
            f = sp.Poly(f, c)
            if f.degree() == 0: continue
            fc = f.as_expr(); matched = False
            for basis_expr, idx in [(c,1),(1-c,2),(1+c,3)]:
                ratio = sp.simplify(fc/basis_expr)
                if ratio.is_number:
                    r = sp.Rational(ratio); r = -r if r<0 else r
                    for prime, ee in sp.factorint(int(r.p)).items():
                        if prime != 2: return None
                        vec[0] += sgn*e*ee
                    for prime, ee in sp.factorint(int(r.q)).items():
                        if prime != 2: return None
                        vec[0] -= sgn*e*ee
                    vec[idx] += sgn*e; matched = True; break
            if not matched: return None
    return tuple(vec)

vals = {}
for a in [-3,-2,-1,0,1,2,3]:
    for b in [-2,-1,0,1,2]:
        for dd in [-2,-1,0,1,2]:
            for e in [-2,-1,0,1,2]:
                base = sp.cancel(2**a * c**b * (1-c)**dd * (1+c)**e)
                for t in (base, -base):
                    if t in (0,1,-1) or t in vals: continue
                    vt = factor_vec(t)
                    if vt is None: continue
                    v1t = factor_vec(1-t)
                    if v1t is None: continue
                    vals[t] = (vt, v1t)
vals[sp.Integer(1)] = ((0,0,0,0),(1,0,0,0))
vals[sp.Integer(-1)] = ((0,0,0,0),(1,0,0,0))
print("S-unit args (deg<=2 grid):", len(vals))
import pickle
pickle.dump({str(k): v for k,v in vals.items()}, open('sunit_args_fixed.pkl','wb'))
# distinct letter pairs for Li2/Li3/Li4 (drop [2]: project vectors to c,p,q coords 1..3)
def proj(v): return v[1:]
seen = {}
for t,(vt,v1t) in vals.items():
    seen.setdefault((proj(v1t), proj(vt)), str(t))   # ([1-g],[g])
print("distinct ([1-g],[g]) pairs:", len(seen))
pickle.dump(seen, open('li_arg_pairs.pkl','wb'))
print(sorted(seen.items())[:10])
