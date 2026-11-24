import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
from dod3_symbol import *
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+mp.sqrt(3)); z2,z3c = mp.pi**2/6, mp.zeta(3)
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
val = pickle.load(open('selftest_val.pkl','rb'))['val']
rel = mp.pslq([val, L3c, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-600'), maxcoeff=10**14, maxsteps=20000)
print("PSLQ [val, L3chi, z3, z2l2, z2l3, z2lu]:", rel)
if rel:
    resid = abs(rel[0]*val + rel[1]*L3c + rel[2]*z3c + rel[3]*z2*l2 + rel[4]*z2*l3 + rel[5]*z2*lu)
    print("resid:", mp.nstr(resid,3))
else:
    # look at the magnitude: maybe coefficient of L3c isn't -1; scan c*L3c subtract
    for c in range(-6,7):
        v2 = val + c*L3c
        rel2 = mp.pslq([v2, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-600'), maxcoeff=10**14, maxsteps=20000)
        print("c=%d:"%c, rel2 if rel2 else "none")
