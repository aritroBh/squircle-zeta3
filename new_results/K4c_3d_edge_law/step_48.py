
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq
from fractions import Fraction

JA = json.load(open("modelA_h96.json")); SA = [mpf(s) for s in JA["Sdk"]]
S0,S1,S2,S3,S4 = SA
l2,l3 = mp.log(2), mp.log(3)
pi2 = mp.pi**2; z3 = mp.zeta(3)
sq5 = mp.sqrt(5); sq2 = mp.sqrt(2)
logs2 = [mpf(1), l2, l3, l2**2, l2*l3, l3**2, pi2]
nlogs2 = ["1","l2","l3","l2^2","l2l3","l3^2","pi2"]
# S2 on {1,sq5} x logs2
B = [a*b for a in [mpf(1), sq5] for b in logs2]
N = [f"{an}*{bn}" for an in ["1","sq5"] for bn in nlogs2]
for tol in [mpf('1e-44')]:
    rel, res = pslq([S2]+B, tol=tol)
    print("S2 {1,sq5}x{logs2}: rel =", rel)
    print("   resid:", mp.nstr(res,3) if res else None)
    if rel:
        # independent residual eval
        rr = rel[0]*S2 + sum(rel[i+1]*B[i] for i in range(len(B)))
        print("   direct residual:", mp.nstr(rr,3))
# Also S2 on {1} x logs2 and {1,sq2,sq5,sq10} x small
B2 = logs2[:]
rel, res = pslq([S2]+B2, tol=mpf('1e-44'))
print("S2 on {1}x{logs2}:", rel, "resid:", mp.nstr(res,3) if res else None)
# S3: basis {1,sq5} x {1,l2,l3,l2^2,l2l3,l3^2,pi2,z3}
logs3 = [mpf(1), l2, l3, l2**2, l2*l3, l3**2, pi2, z3]
nlogs3 = ["1","l2","l3","l2^2","l2l3","l3^2","pi2","z3"]
B3 = [a*b for a in [mpf(1), sq5] for b in logs3]
rel, res = pslq([S3]+B3, tol=mpf('1e-44'))
print("S3 {1,sq5}x{logs3}:", rel, "resid:", mp.nstr(res,3) if res else None)
if rel:
    rr = rel[0]*S3 + sum(rel[i+1]*B3[i] for i in range(len(B3)))
    print("   direct residual:", mp.nstr(rr,3))
    a0=rel[0]
    combos=[(f"{['1','sq5'][i//8]}*{nlogs3[i%8]}", Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(B3)) if rel[i+1]!=0]
    print("   S3 =", " + ".join(f"({c})*{n}" for n,c in combos))
# S3 on {1} x logs3
rel, res = pslq([S3]+logs3, tol=mpf('1e-44'))
print("S3 on {1}x{logs3}:", rel, "resid:", mp.nstr(res,3) if res else None)
