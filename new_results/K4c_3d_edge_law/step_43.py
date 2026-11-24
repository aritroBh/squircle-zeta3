
with open('pslq_lib.py') as f: txt = f.read()
txt = txt.replace("for j in range(min(i, n-2), -1, -1):", "for j in range(min(i-1, n-2), -1, -1):")
with open('pslq_lib.py','w') as f: f.write(txt)
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
import importlib, pslq_lib
importlib.reload(pslq_lib)
from pslq_lib import pslq
q0 = mpf(1)/2 + 2*mp.log(2)
rel, res = pslq([q0, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-45'))
print("synthetic q0 (1,ln2,ln3):", rel, "resid:", mp.nstr(res,3))
rel, res = pslq([mp.log(4), mp.log(2), mpf(1)], tol=mpf('1e-45'))
print("ln4 test:", rel)
import json
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
S = [mpf(s_) for s_ in J["Sdk"]]
q = S[1]/S[0]
rel, res = pslq([q, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-45'))
print("cubo q:", rel, "resid:", mp.nstr(res,3))
for name, fn in [("A","modelA_h96.json"), ("B","modelB_h48.json")]:
    J2 = json.load(open(fn)); S2 = [mpf(s_) for s_ in J2["Sdk"]]
    q2 = S2[1]/S2[0]
    for bname, basis in [("(1,ln2,ln3)", [mpf(1), mp.log(2), mp.log(3)]),
                         ("(1,ln2,ln3,ln5)", [mpf(1), mp.log(2), mp.log(3), mp.log(5)])]:
        rel, res = pslq([q2]+basis, tol=mpf('1e-40'))
        print(f"model {name} q on {bname}:", rel, "resid:", mp.nstr(res,3) if res else None)
