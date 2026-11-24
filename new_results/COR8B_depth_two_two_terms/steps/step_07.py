
# step_07.py -- expanded S-unit enumeration (signed monomials) + delta22 database
import sys, json, pickle, time
sys.path.insert(0,'.')
from cor8b_sym import *
import sympy as sp
c = sp.symbols('c')

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
exps = [-3,-2,-1,0,1,2,3]
exps2 = [-2,-1,0,1,2]
t0=time.time()
for a in exps2:
    for b in exps:
        for dd in exps:
            for e in exps:
                base = sp.cancel(2**a * c**b * (1-c)**dd * (1+c)**e)
                for t in (base, -base):
                    if t in (0,1,-1): continue
                    if t in vals: continue
                    vt = factor_vec(t)
                    if vt is None: continue
                    v1t = factor_vec(1-t)
                    if v1t is None: continue
                    vals[t] = (vt, v1t)
print("valid S-unit args (signed):", len(vals), "time", time.time()-t0)
# add 1 and -1 explicitly as arguments (letters [1]=0, [1-t]=[2] or [0])
vals[sp.Integer(1)] = ((0,0,0,0), (1,0,0,0))
vals[sp.Integer(-1)] = ((0,0,0,0), (1,0,0,0))  # [1-(-1)]=[2]
pairs = []
vl = list(vals.items())
t0=time.time()
for x0,(vx,vux) in vl:
    for y0,(vy,vyv) in vl:
        xy = sp.cancel(x0*y0)
        if xy == 1:
            vz, vw = (0,0,0,0), None
            # 1-xy = 0 -> W zero vector allowed
            vw = (0,0,0,0)
            pairs.append((x0,y0,vx,vux,vy,vyv,vz,vw)); continue
        vw = factor_vec(1-xy)
        if vw is None: continue
        vz = factor_vec(xy)
        pairs.append((x0,y0,vx,vux,vy,vyv,vz,vw))
print("valid pairs:", len(pairs), "time", time.time()-t0)
S22 = S_Li22(); S31 = S_Li31(); S13 = S_Li13()
t0=time.time()
data = []
for (x0,y0,vx,vux,vy,vyv,vz,vw) in pairs:
    imgs = {'X':vx,'U':vux,'Y':vy,'V':vyv,'W':vw}
    d22 = delta22(specialize_pair(S22, imgs))
    d31 = delta22(specialize_pair(S31, imgs))
    d13 = delta22(specialize_pair(S13, imgs))
    if d22 or d31 or d13:
        data.append((str(x0),str(y0),d22,d31,d13))
print("pairs with nonzero delta:", len(data), "time", time.time()-t0)
pickle.dump(data, open('pairs_delta_all_v2.pkl','wb'))
from collections import defaultdict
byvec = defaultdict(int)
for (_,_,d22,d31,d13) in data:
    for dd in (d22,d31,d13):
        byvec[tuple(sorted(dd.items()))]+=1
print("distinct delta vectors:", len(byvec))
