
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq

JA = json.load(open("modelA_h96.json")); SA = [mpf(s) for s in JA["Sdk"]]
S0,S1,S2,S3,S4 = SA
L = mpf(-6)/5*mp.log(2) + mpf(13)/5*mp.log(3)
e1 = L
e2 = L**2/2 - L
e3 = L**3/6 - L**2 + L
e4 = L**4/24 - L**3/2 + 3*L**2/2 - L
t0 = S0
t1 = S1 - e1*t0
t2 = S2 - e1*t1 - e2*t0
t3 = S3 - e1*t2 - e2*t1 - e3*t0
t4 = S4 - e1*t3 - e2*t2 - e3*t1 - e4*t0
print("t1/S0 =", mp.nstr(t1/S0,30), " (expect 13/15 =", mp.nstr(mpf(13)/15,30),")")
print("t2 =", mp.nstr(t2,50))
print("t3 =", mp.nstr(t3,50))
print("t4 =", mp.nstr(t4,50))
sq = {r: mp.sqrt(r) for r in [2,3,5,6]}
pi2 = mp.pi**2; z3 = mp.zeta(3)
# basis: {1,sq2,sq3,sq5,sq6} x {1,pi2,z3}
basis = []
names = []
for rn, rv in [("1",mpf(1)), ("sq2",sq[2]), ("sq3",sq[3]), ("sq5",sq[5]), ("sq6",sq[6])]:
    for cn, cv in [("1",mpf(1)), ("pi2",pi2), ("z3",z3)]:
        basis.append(rv*cv); names.append(f"{rn}*{cn}")
for tn, tv in [("t2",t2), ("t3",t3)]:
    rel, res = pslq([tv]+basis, tol=mpf('1e-40'))
    if rel:
        print(f"{tn}: rel={rel} resid={mp.nstr(res,3)}")
        # print as combination
        a0 = rel[0]
        from fractions import Fraction
        combos = [(names[i], Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(basis)) if rel[i+1]!=0]
        print("   ", tn, "=", " + ".join(f"({c})*{n}" for n,c in combos))
    else:
        print(f"{tn}: no relation at 1e-40")
