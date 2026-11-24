# step_85: find Kummer-Spence sign pattern with delta=0 (fixed wedge), then test
from fractions import Fraction as Fr
import itertools
def v(**kw): return {k:Fr(c) for k,c in kw.items() if c}
pts = {
 'x':      (v(x=1),          v(X=1)),
 'y':      (v(y=1),          v(Y=1)),
 'y/x':    (v(y=1,x=-1),     v(D=1,x=-1)),
 '(1-x)/(1-y)':  (v(X=1,Y=-1),  v(D=1,Y=-1)),
 '(1-x)/(1-xy)': (v(X=1,W=-1),  v(x=1,y=1,W=-1)),
 '(1-y)/(1-xy)': (v(Y=1,W=-1),  v(x=1,y=1,W=-1)),
 'x(1-y)/(1-xy)': (v(x=1,Y=1,W=-1), v(X=1,W=-1)),
 'y(1-x)/(1-xy)': (v(y=1,X=1,W=-1), v(Y=1,W=-1)),
}
names = list(pts.keys())
def wedge(u,w):
    out={}
    for a,ca in u.items():
        for b,cb in w.items():
            if a==b: continue
            if a<b: key, s = (a,b), 1
            else: key, s = (b,a), -1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_of_pattern(eps):
    out={}
    for nm,e in zip(names,eps):
        if e==0: continue
        vz, vomz = pts[nm]
        w = wedge(vz,vomz)
        for (a,b),cab in w.items():
            for l3,c3 in vz.items():
                key=((a,b),l3)
                out[key]=out.get(key,Fr(0))+e*cab*c3
    return {k:c for k,c in out.items() if c}
sols = []
for eps in itertools.product([1,-1], repeat=8):
    if not delta_of_pattern(eps):
        sols.append(eps)
print("sign patterns with delta = 0:", sols)
import mpmath as mp
mp.mp.dps = 100
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
for eps in sols:
    for xx, yy in [('0.3','0.6'), ('0.2','0.7'), ('0.45','0.85')]:
        x, y = mp.mpf(xx), mp.mpf(yy)
        args = [x, y, y/x, (1-x)/(1-y), (1-x)/(1-x*y), (1-y)/(1-x*y), x*(1-y)/(1-x*y), y*(1-x)/(1-x*y)]
        val = mp.fsum([e*L3(z) for e,z in zip(eps, args)])
        print(eps, f"({xx},{yy}):", mp.nstr(val, 8))
