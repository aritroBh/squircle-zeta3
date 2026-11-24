
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

t0,t1,t2,t3 = get_tr("modelA_h96.json")
pi2=mp.pi**2; z3=mp.zeta(3)
sq2,sq3,sq5,sq6 = mp.sqrt(2),mp.sqrt(3),mp.sqrt(5),mp.sqrt(6)
basis15 = [mpf(1),sq2,sq3,sq5,sq6, pi2,sq2*pi2,sq3*pi2,sq5*pi2,sq6*pi2, z3,sq2*z3,sq3*z3,sq5*z3,sq6*z3]
n15 = ["1","sq2","sq3","sq5","sq6","pi2","sq2*pi2","sq3*pi2","sq5*pi2","sq6*pi2","z3","sq2*z3","sq3*z3","sq5*z3","sq6*z3"]
# candidate relations from before
rel_t2 = [168, -461, -68, 255, -98, -211, -188, -30, 165, -79, -145, 112, 155, -186, 28, -367]
resid = rel_t2[0]*t2 + sum(rel_t2[i+1]*basis15[i] for i in range(15))
print("t2 candidate residual (h96 data):", mp.nstr(resid,5))
rel_t3 = [131, 18, 180, -160, 6, 67, -152, 89, 112, 27, 157, -246, -230, 149, 30, 82]
resid = rel_t3[0]*t3 + sum(rel_t3[i+1]*basis15[i] for i in range(15))
print("t3 candidate residual (h96 data):", mp.nstr(resid,5))
# same on h64 and h48
t0_64,t1_64,t2_64,t3_64 = get_tr("modelA_h64.json")
resid = rel_t2[0]*t2_64 + sum(rel_t2[i+1]*basis15[i] for i in range(15))
print("t2 candidate residual (h64 data):", mp.nstr(resid,5))
resid = rel_t3[0]*t3_64 + sum(rel_t3[i+1]*basis15[i] for i in range(15))
print("t3 candidate residual (h64 data):", mp.nstr(resid,5))
# PSLQ on h64 data (independent): does it find the same?
rel2, res2 = pslq([t2_64]+basis15, tol=mpf('1e-30'))
print("t2 PSLQ on h64:", rel2, mp.nstr(res2,3) if res2 else None)
rel3, res3 = pslq([t3_64]+basis15, tol=mpf('1e-30'))
print("t3 PSLQ on h64:", rel3, mp.nstr(res3,3) if res3 else None)
