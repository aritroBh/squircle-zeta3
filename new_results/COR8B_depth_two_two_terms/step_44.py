
script = r'''
# step_08.py -- randomized robustness test of the 4 obstruction constraints (large exponents)
import sys, random, time
sys.path.insert(0,'.')
from cor8b_sym import *
import sympy as sp
c = sp.symbols('c')
random.seed(12345)

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

S22 = S_Li22(); S31 = S_Li31(); S13 = S_Li13()
USE = [(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)]
uidx = {k:i for i,k in enumerate(USE)}
CONS = [((1,3),(3,4)), ((1,4),(4,5))]  # not used; we use 4 functionals below
def vec_of(dd):
    v = [0]*7
    for k,val in dd.items():
        kk=(min(k),max(k)); sgn=1 if k==kk else -1
        if kk in uidx: v[uidx[kk]] += sgn*val
    return v
def cons_ok(v):
    # f1 mp^cp+mq^cp = v0+v2 ; f2 mp^cq+mq^cq = v1+v3 ; f3 mp^cp+cp^pq = v0+v5 ; f4 mp^cq+cq^pq = v1+v6
    return (v[0]+v[2]==0) and (v[1]+v[3]==0) and (v[0]+v[5]==0) and (v[1]+v[6]==0)

valid_args = []
t0=time.time()
ntries=0
while len(valid_args) < 400 and time.time()-t0 < 500:
    ntries += 1
    a = random.randint(-3,3); b = random.randint(-6,6); dd = random.randint(-6,6); e = random.randint(-6,6)
    base = sp.cancel(2**a * c**b * (1-c)**dd * (1+c)**e)
    t = base if random.random()<0.5 else -base
    if t in (0,1,-1): continue
    vt = factor_vec(t)
    if vt is None: continue
    v1t = factor_vec(1-t)
    if v1t is None: continue
    valid_args.append((t, vt, v1t))
print("valid args from sampling:", len(valid_args), "tries", ntries)

tested = 0; bad = 0; examples = []
t0=time.time()
for _ in range(3000):
    x0, vx, vux = random.choice(valid_args)
    y0, vy, vyv = random.choice(valid_args)
    xy = sp.cancel(x0*y0)
    if xy == 1: vz, vw = (0,0,0,0),(0,0,0,0)
    else:
        vw = factor_vec(1-xy)
        if vw is None: continue
        vz = factor_vec(xy)
    imgs = {'X':vx,'U':vux,'Y':vy,'V':vyv,'W':vw}
    for S in (S22,S31,S13):
        dd = delta22(specialize_pair(S, imgs))
        v = vec_of(dd)
        tested += 1
        if not cons_ok(v):
            bad += 1
            examples.append((str(x0),str(y0),v))
print("tested:", tested, "violations:", bad, "time", time.time()-t0)
for ex in examples[:10]: print("  VIOLATION:", ex)
open('step08_result.txt','w').write(f"tested {tested} violations {bad}\n" + "\n".join(map(str,examples[:20])))
'''
open('steps/step_08.py','w').write(script)
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','15','<project>/.venv/bin/python','steps/step_08.py'],capture_output=True,text=True,timeout=1400)
print(r.stdout[-3000:]); print(r.stderr[-2000:]); print("elapsed", time.time()-t)
open('progress.log','a').write("step_08 randomized obstruction robustness test\n")
