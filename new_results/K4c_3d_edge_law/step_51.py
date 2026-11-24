
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 65
from pslq_lib import pslq
JB = json.load(open("modelB_h48.json")); SB = [mpf(s) for s in JB["Sdk"]]
qB = SB[1]/SB[0]
sq3 = mp.sqrt(3)
Cl2 = mp.clsin(mp.pi/3, 2) if hasattr(mp,'clsin') else None
# Clausen Cl2(pi/3) via mpmath: mp.clsin(s, z) gives sum sin(k s)/k^z? clsin(theta, q) = Im Li_q(e^{i theta}) for q=2 -> Cl2(theta)
print("Cl2(pi/3) =", mp.nstr(Cl2,30) if Cl2 else "n/a")
G = mp.catalan
bases = [
 ("1,l2,l3,l7,Cl2", [mpf(1), mp.log(2), mp.log(3), mp.log(7), Cl2]),
 ("1,l2,l3,l7,Cl2,G", [mpf(1), mp.log(2), mp.log(3), mp.log(7), Cl2, G]),
 ("1,l2,l3,l7,l13,Cl2", [mpf(1), mp.log(2), mp.log(3), mp.log(7), mp.log(13), Cl2]),
 ("1,l2,l3,l7,l13,pi", [mpf(1), mp.log(2), mp.log(3), mp.log(7), mp.log(13), mp.pi]),
 ("1,l2,l3,l7", [mpf(1), mp.log(2), mp.log(3), mp.log(7)]),
 ("1,l2,l3,l(2+sq3),Cl2", [mpf(1), mp.log(2), mp.log(3), mp.log(2+sq3), Cl2]),
 ("l2,l3,l7,l13,l(2+sq3)", [mp.log(2), mp.log(3), mp.log(7), mp.log(13), mp.log(2+sq3)]),
]
for bn, b in bases:
    if None in b: continue
    rel, res = pslq([qB]+b, tol=mpf('1e-50'))
    sig = ""
    if rel:
        height = max(abs(t) for t in rel)
        sig = "height=%g" % height
    print(f"qB on {bn}: {rel} resid={mp.nstr(res,3) if res else None} {sig}")
