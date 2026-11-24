
script = r'''
# step_06.py -- enumerate S-unit arguments and compute all depth-two delta22's
import sys, json, pickle, time
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr
import sympy as sp

c = sp.symbols('c')

def factor_vec(expr):
    # expr: sympy rational function in c. Return 4-vector over (2,c,1-c,1+c) or None.
    expr = sp.simplify(expr)
    num, den = sp.fraction(sp.cancel(expr))
    vec = [0,0,0,0]
    for poly, sgn in [(num,1),(den,-1)]:
        poly = sp.Poly(poly, c)
        content, factors = sp.factor_list(poly)
        # content is a rational number
        cont = sp.Rational(content)
        if cont < 0: cont = -cont  # [-1]=0
        for pr in sp.factorint(abs(cont.p)).items(): pass
        for prime, e in sp.factorint(int(abs(cont.p*cont.q) if False else cont.p)).items():
            if prime != 2: return None
            vec[0] += sgn*e
        for prime, e in sp.factorint(int(cont.q)).items():
            if prime != 2: return None
            vec[0] -= sgn*e
        for f, e in factors:
            f = sp.Poly(f, c)
            if f.degree() == 0:
                continue
            # normalize: compare with c, 1-c, 1+c up to scalar
            fc = f.as_expr()
            matched = False
            for basis_expr, idx in [(c,1),(1-c,2),(1+c,3)]:
                ratio = sp.simplify(fc/basis_expr)
                if ratio.is_number:
                    r = sp.Rational(ratio)
                    if r < 0: r = -r
                    for prime, ee in sp.factorint(int(r.p)).items():
                        if prime != 2: return None
                        vec[0] += sgn*e*ee
                    for prime, ee in sp.factorint(int(r.q)).items():
                        if prime != 2: return None
                        vec[0] -= sgn*e*ee
                    vec[idx] += sgn*e
                    matched = True
                    break
            if not matched:
                return None
    return tuple(vec)

# monomials 2^a c^b (1-c)^d (1+c)^e
vals = set()
exps = [-2,-1,0,1,2]
for a in exps:
    for b in exps:
        for dd in exps:
            for e in exps:
                t = 2**a * c**b * (1-c)**dd * (1+c)**e
                t = sp.cancel(t)
                if t in (0,1): continue
                vt = factor_vec(t)
                if vt is None: continue
                v1t = factor_vec(1-t)
                if v1t is None: continue
                vals.add((sp.cancel(t), vt, v1t))
print("valid S-unit args:", len(vals))
# pairs
pairs = []
vl = list(vals)
t0=time.time()
for i,(x0,vx,vux) in enumerate(vl):
    for (y0,vy,vyv) in vl:
        xy = sp.cancel(x0*y0)
        vw = factor_vec(1-xy)
        if vw is None: continue
        vz = factor_vec(xy)
        pairs.append((x0,y0,vx,vux,vy,vyv,vz,vw))
print("valid pairs:", len(pairs), "time", time.time()-t0)

# compute delta22 for Li22, Li31, Li13 at each pair
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
pickle.dump(data, open('pairs_delta_all.pkl','wb'))
# quick stats: which components appear
comp = {}
for (_,_,d22,d31,d13) in data:
    for dd in (d22,d31,d13):
        for k in dd: comp[k] = comp.get(k,0)+1
print("component usage:", {WNAME[k[0]]+'^'+WNAME[k[1]]: v for k,v in sorted(comp.items())})
'''
open('steps/step_06.py','w').write(script)
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','15','<project>/.venv/bin/python','steps/step_06.py'],capture_output=True,text=True,timeout=1400)
print(r.stdout[-4000:]); print(r.stderr[-2000:]); print("elapsed", time.time()-t)
open('progress.log','a').write("step_06 S-unit pairs + delta22 database\n")
