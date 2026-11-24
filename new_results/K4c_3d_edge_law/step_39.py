
import re
with open('pslq_lib.py') as f: txt = f.read()
txt = txt.replace("""            rel = [B[j][imin] for j in range(n)]
            g = rel[0]
            for t in rel[1:]:
                from math import gcd
                g = gcd(int(g), int(t))
            if g != 0:
                rel = [t//g for t in rel]
            return [int(t) for t in rel], abs(y[imin])""",
"""            rel = [int(mp.nint(B[j][imin])) for j in range(n)]
            from math import gcd
            g = 0
            for t in rel:
                g = gcd(g, abs(t))
            if g > 1:
                rel = [t//g for t in rel]
            return rel, abs(y[imin])""")
with open('pslq_lib.py','w') as f: f.write(txt)
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 50
import importlib, pslq_lib
importlib.reload(pslq_lib)
from pslq_lib import pslq
rel, res = pslq([mp.log(4), mp.log(2), mpf(1)], tol=mpf('1e-40'))
print("test1 (ln4,ln2,1):", rel, "resid:", mp.nstr(res,3))
rel, res = pslq([mp.zeta(3), mp.pi**3, mpf(1)], tol=mpf('1e-38'))
print("test2 (zeta3,pi^3,1):", rel, "resid:", mp.nstr(res,3) if res else None)
import json
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
S = [mpf(s_) for s_ in J["Sdk"]]
q = S[1]/S[0]
rel, res = pslq([q, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-40'))
print("test3 cubo q=1.8862944 on (1,ln2,ln3):", rel, "resid:", mp.nstr(res,3))
for name, fn in [("A","modelA_h96.json"), ("B","modelB_h48.json")]:
    J2 = json.load(open(fn)); S2 = [mpf(s_) for s_ in J2["Sdk"]]
    q2 = S2[1]/S2[0]
    for bname, basis in [("(1,ln2,ln3)", [mpf(1), mp.log(2), mp.log(3)]),
                         ("(1,ln2,ln3,ln5)", [mpf(1), mp.log(2), mp.log(3), mp.log(5)]),
                         ("(1,ln2,ln3,ln5,ln7)", [mpf(1), mp.log(2), mp.log(3), mp.log(5), mp.log(7)])]:
        rel, res = pslq([q2]+basis, tol=mpf('1e-38'))
        print(f"model {name} q on {bname}:", rel, "resid:", mp.nstr(res,3) if res else None)
