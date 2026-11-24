
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq
from fractions import Fraction

JA = json.load(open("modelA_h96.json")); SA = [mpf(s) for s in JA["Sdk"]]
S0,S1,S2,S3,S4 = SA
L = mpf(-6)/5*mp.log(2) + mpf(13)/5*mp.log(3)
e1=L; e2=L**2/2-L; e3=L**3/6-L**2+L; e4=L**4/24-L**3/2+3*L**2/2-L
t0=S0; t1=S1-e1*t0; t2=S2-e1*t1-e2*t0; t3=S3-e1*t2-e2*t1-e3*t0; t4=S4-e1*t3-e2*t2-e3*t1-e4*t0
r2, r3 = t2/S0, t3/S0
print("r2 =", mp.nstr(r2,50))
print("r3 =", mp.nstr(r3,50))
pi2=mp.pi**2; z3=mp.zeta(3)
sq2,sq3,sq5,sq6 = mp.sqrt(2),mp.sqrt(3),mp.sqrt(5),mp.sqrt(6)

def try_basis(vals, basis, names, tol):
    rel, res = pslq([vals]+basis, tol=tol)
    if rel:
        a0 = rel[0]
        combos = [(names[i], Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(basis)) if rel[i+1]!=0]
        print(f"   rel resid={mp.nstr(res,3)}: {tn} =", " + ".join(f"({c})*{n}" for n,c in combos))
    else:
        print(f"   none at {mp.nstr(tol,2)}")

b6 = [mpf(1), sq2, pi2, sq2*pi2, z3, sq2*z3]
n6 = ["1","sq2","pi2","sq2*pi2","z3","sq2*z3"]
b12 = []
n12 = []
for rn, rv in [("1",mpf(1)),("sq2",sq2),("sq3",sq3),("sq6",sq6)]:
    for cn, cv in [("1",mpf(1)),("pi2",pi2),("z3",z3)]:
        b12.append(rv*cv); n12.append(f"{rn}*{cn}")
for tn, tv in [("r2",r2),("r3",r3)]:
    print(f"{tn} on {{1,sq2}}x{{1,pi2,z3}}:")
    try_basis(tv, b6, n6, mpf('1e-42'))
    print(f"{tn} on {{1,sq2,sq3,sq6}}x{{1,pi2,z3}}:")
    try_basis(tv, b12, n12, mpf('1e-40'))
