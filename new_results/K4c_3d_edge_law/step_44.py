
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq

JA = json.load(open("modelA_h96.json")); SA = [mpf(s) for s in JA["Sdk"]]
JB = json.load(open("modelB_h48.json")); SB = [mpf(s) for s in JB["Sdk"]]
qA, qB = SA[1]/SA[0], SB[1]/SB[0]
print("qA =", mp.nstr(qA,55))
print("qB =", mp.nstr(qB,55))
# verify model A relation at dps 70
rel, res = pslq([qA, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-50'))
print("qA on (1,ln2,ln3):", rel, "resid:", mp.nstr(res,3))
# model B extended bases
sq3 = mp.sqrt(3)
cands = {
 "ln2": mp.log(2), "ln3": mp.log(3), "ln13": mp.log(13),
 "ln(2+sq3)": mp.log(2+sq3), "ln(4sq3-3)": mp.log(4*sq3-3),
 "ln(2sq3-3)": mp.log(2*sq3-3), "ln7": mp.log(7), "ln(3+2sq3)": mp.log(3+2*sq3),
 "ln(433)": mp.log(433),
}
bases = [
  ("1,l2,l3,l13", [mpf(1), cands["ln2"], cands["ln3"], cands["ln13"]]),
  ("1,l2,l3,l(2+sq3)", [mpf(1), cands["ln2"], cands["ln3"], cands["ln(2+sq3)"]]),
  ("1,l2,l3,l(2sq3-3)", [mpf(1), cands["ln2"], cands["ln3"], cands["ln(2sq3-3)"]]),
  ("1,l2,l3,l13,l(2+sq3)", [mpf(1), cands["ln2"], cands["ln3"], cands["ln13"], cands["ln(2+sq3)"]]),
  ("1,l2,l3,l13,l(2sq3-3),l(4sq3-3)", [mpf(1), cands["ln2"], cands["ln3"], cands["ln13"], cands["ln(2sq3-3)"], cands["ln(4sq3-3)"]]),
]
for bn, b in bases:
    rel, res = pslq([qB]+b, tol=mpf('1e-45'))
    print(f"qB on {bn}:", rel, "resid:", mp.nstr(res,3) if res else None)
