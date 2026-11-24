
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq
from fractions import Fraction

def get_tr(fn):
    J = json.load(open(fn)); S = [mpf(s) for s in J["Sdk"]]
    S0,S1,S2,S3,S4 = S
    L = mpf(-6)/5*mp.log(2) + mpf(13)/5*mp.log(3)
    e1=L; e2=L**2/2-L; e3=L**3/6-L**2+L
    t0=S0; t1=S1-e1*t0; t2=S2-e1*t1-e2*t0; t3=S3-e1*t2-e2*t1-e3*t0
    return t0,t1,t2,t3

pi2=mp.pi**2; z3=mp.zeta(3)
sq2,sq3,sq5,sq6,sq10 = mp.sqrt(2),mp.sqrt(3),mp.sqrt(5),mp.sqrt(6),mp.sqrt(10)
b12 = []
n12 = []
for rn, rv in [("1",mpf(1)),("sq2",sq2),("sq3",sq3),("sq6",sq6)]:
    for cn, cv in [("1",mpf(1)),("pi2",pi2),("z3",z3)]:
        b12.append(rv*cv); n12.append(f"{rn}*{cn}")

for fn, tag in [("modelA_h96.json","h96"), ("modelA_h64.json","h64"), ("modelA_h48.json","h48")]:
    t0,t1,t2,t3 = get_tr(fn)
    r2, r3 = t2/t0, t3/t0
    rel2, res2 = pslq([r2]+b12, tol=mpf('1e-38'))
    rel3, res3 = pslq([r3]+b12, tol=mpf('1e-38'))
    print(f"{tag}: r2 rel={rel2} resid={mp.nstr(res2,3) if res2 else None}")
    print(f"{tag}: r3 rel={rel3} resid={mp.nstr(res3,3) if res3 else None}")
# also {1,sq5} x {1,pi2,z3} for r2,r3 (6-dim)
b6 = [mpf(1), sq5, pi2, sq5*pi2, z3, sq5*z3]
for fn, tag in [("modelA_h96.json","h96"), ("modelA_h64.json","h64")]:
    t0,t1,t2,t3 = get_tr(fn)
    r2, r3 = t2/t0, t3/t0
    rel2, res2 = pslq([r2]+b6, tol=mpf('1e-42'))
    rel3, res3 = pslq([r3]+b6, tol=mpf('1e-42'))
    print(f"{tag} {{1,sq5}}x: r2 rel={rel2} resid={mp.nstr(res2,3) if res2 else None}")
    print(f"{tag} {{1,sq5}}x: r3 rel={rel3} resid={mp.nstr(res3,3) if res3 else None}")
