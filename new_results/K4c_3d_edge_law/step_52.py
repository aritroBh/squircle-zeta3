
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 65
from pslq_lib import pslq
JB = json.load(open("modelB_h48.json")); SB = [mpf(s) for s in JB["Sdk"]]
qB = SB[1]/SB[0]
sq3 = mp.sqrt(3); sq2 = mp.sqrt(2)
consts = {
 "l2": mp.log(2), "l3": mp.log(3), "l7": mp.log(7), "l13": mp.log(13),
 "l(2+s3)": mp.log(2+sq3), "l(1+s2)": mp.log(1+sq2), "l(s2+s3)": mp.log(sq2+sq3),
 "lh_sq": mp.log(1+4/sq3), "lh_hex": mp.log(sq2*(sq3+2)/sq3),
 "l(4s3-3)": mp.log(4*sq3-3), "l(2s3-3)": mp.log(2*sq3-3),
}
tests = [
 ["1","l2","l3","l(2+s3)","lh_sq","lh_hex"],
 ["1","l2","l3","l7","lh_sq","lh_hex"],
 ["1","l2","l3","l13","l(4s3-3)","l(2s3-3)"],
 ["1","l2","l3","lh_sq"],
 ["1","l2","l3","lh_hex"],
 ["1","l2","l3","l(2+s3)","l(1+s2)","l(s2+s3)"],
 ["l2","l3","l7","lh_sq","lh_hex"],
 ["1","l2","l3","l7","l13","l(2+s3)"],
]
bmap = {"1": mpf(1)}
for t in tests:
    b = [bmap.get(x, consts.get(x)) for x in t]
    rel, res = pslq([qB]+b, tol=mpf('1e-48'))
    h = max(abs(x) for x in rel) if rel else None
    print(f"qB on {t}: height={h} resid={mp.nstr(res,3) if res else None}")
    if rel and h < 10000:
        print("   rel:", rel)
# Also S1 directly on {1, sq3} x {1, l2, l3, l7, l13}
S1 = SB[1]
B = []
N = []
for an, av in [("1",mpf(1)), ("sq3",sq3)]:
    for cn, cv in [("1",mpf(1)), ("l2",mp.log(2)), ("l3",mp.log(3)), ("l7",mp.log(7)), ("l13",mp.log(13))]:
        B.append(av*cv); N.append(f"{an}*{cn}")
rel, res = pslq([S1]+B, tol=mpf('1e-48'))
print("S1(B) on {1,sq3}x{1,l2,l3,l7,l13}:", rel, "resid:", mp.nstr(res,3) if res else None)
