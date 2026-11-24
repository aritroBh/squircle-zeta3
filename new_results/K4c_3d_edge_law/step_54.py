
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq
from fractions import Fraction
JA = json.load(open("modelA_h96.json")); SA = [mpf(s) for s in JA["Sdk"]]
S0,S1,S2,S3,S4 = SA
sq5 = mp.sqrt(5)
l2,l3 = mp.log(2), mp.log(3)
pi2, pi3, z3 = mp.pi**2, mp.pi**3, mp.zeta(3)

# S2/sqrt5 on rational log basis
b2 = [mpf(1), l2, l3, l2**2, l2*l3, l3**2, pi2]
n2 = ["1","l2","l3","l2^2","l2l3","l3^2","pi2"]
v = S2/sq5
rel, res = pslq([v]+b2, tol=mpf('1e-46'))
print("S2/sq5 on {1,l2,l3,l2^2,l2l3,l3^2,pi2}:", rel, "resid:", mp.nstr(res,3) if res else None)
if rel:
    a0=rel[0]
    combos=[(n2[i], Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(b2)) if rel[i+1]!=0]
    print("   S2/sq5 =", " + ".join(f"({c})*{n}" for n,c in combos))
# S3/sqrt5 on rational log basis (weight <=3)
b3 = [mpf(1), l2, l3, l2**2, l2*l3, l3**2, pi2, l2**3, l2**2*l3, l2*l3**2, l3**3, pi3, pi2*l2, pi2*l3, z3]
n3 = ["1","l2","l3","l2^2","l2l3","l3^2","pi2","l2^3","l2^2l3","l2l3^2","l3^3","pi3","pi2l2","pi2l3","z3"]
v = S3/sq5
rel, res = pslq([v]+b3, tol=mpf('1e-42'))
print("S3/sq5 on 15-dim weight<=3:", rel, "resid:", mp.nstr(res,3) if res else None)
if rel:
    a0=rel[0]
    combos=[(n3[i], Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(b3)) if rel[i+1]!=0]
    print("   S3/sq5 =", " + ".join(f"({c})*{n}" for n,c in combos))
# S1/sqrt5 as sanity (should reproduce q_A * (2/3))
b1 = [mpf(1), l2, l3]
v = S1/sq5
rel, res = pslq([v]+b1, tol=mpf('1e-50'))
print("S1/sq5 on {1,l2,l3}:", rel, "resid:", mp.nstr(res,3) if res else None)
