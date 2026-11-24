
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 65
from pslq_lib import pslq
JB = json.load(open("modelB_h48.json")); SB = [mpf(s) for s in JB["Sdk"]]
qB = SB[1]/SB[0]
sq3 = mp.sqrt(3); sq2 = mp.sqrt(2)
l2, l3, l7, l13 = mp.log(2), mp.log(3), mp.log(7), mp.log(13)
lA = mp.log(2+sq3); lB = mp.log(4+sq3); lC = mp.log(1+sq2); lD = mp.log(sq2+sq3)
tests = [
 ("1,l2,l3,l7,l13,lA,lB", [mpf(1), l2, l3, l7, l13, lA, lB]),
 ("1,l2,l3,lA,lB", [mpf(1), l2, l3, lA, lB]),
 ("1,l2,l3,l7,lA,lB", [mpf(1), l2, l3, l7, lA, lB]),
 ("1,l2,l3,l13,lA,lB", [mpf(1), l2, l3, l13, lA, lB]),
 ("1,l2,l3,l7,l13,lA,lB,lC,lD", [mpf(1), l2, l3, l7, l13, lA, lB, lC, lD]),
]
for bn, b in tests:
    rel, res = pslq([qB]+b, tol=mpf('1e-50'))
    h = max(abs(x) for x in rel) if rel else None
    print(f"qB on {bn}: height={h} resid={mp.nstr(res,3) if res else None}")
    if rel and h < 10**6:
        print("   rel:", rel)
