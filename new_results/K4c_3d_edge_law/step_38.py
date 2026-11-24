
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 50
from pslq_lib import pslq
x = [mp.log(4), mp.log(2), mpf(1)]
rel, res = pslq(x, tol=mpf('1e-40'))
print("test1 (ln4,ln2,1):", rel, "resid:", mp.nstr(res,3))
rel, res = pslq([mp.zeta(3), mp.pi**3, mpf(1)], tol=mpf('1e-38'))
print("test2 (zeta3,pi^3,1):", rel, "resid:", mp.nstr(res,3) if res else None)
import json
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
S = [mpf(s_) for s_ in J["Sdk"]]
q = S[1]/S[0]
rel, res = pslq([q, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-40'))
print("test3 cubo q=1.886... on (1,ln2,ln3):", rel, "resid:", mp.nstr(res,3))
# model A and B q values
for name, fn in [("A","modelA_h96.json"), ("B","modelB_h48.json")]:
    J2 = json.load(open(fn)); S2 = [mpf(s_) for s_ in J2["Sdk"]]
    q2 = S2[1]/S2[0]
    for bname, basis in [("(1,ln2,ln3)", [mpf(1), mp.log(2), mp.log(3)]),
                         ("(1,ln2,ln3,ln5)", [mpf(1), mp.log(2), mp.log(3), mp.log(5)]),
                         ("(1,ln2,ln3,ln5,ln7)", [mpf(1), mp.log(2), mp.log(3), mp.log(5), mp.log(7)])]:
        rel, res = pslq([q2]+basis, tol=mpf('1e-38'))
        print(f"model {name} q on {bname}:", rel, "resid:", mp.nstr(res,3) if res else None)
