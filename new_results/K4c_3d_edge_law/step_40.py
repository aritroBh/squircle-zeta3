
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
import importlib, pslq_lib
importlib.reload(pslq_lib)
from pslq_lib import pslq

# synthetic: q0 = 1/2 + 2 ln2
q0 = mpf(1)/2 + 2*mp.log(2)
for tol in [mpf('1e-30'), mpf('1e-45')]:
    rel, res = pslq([q0, mpf(1), mp.log(2), mp.log(3)], tol=tol)
    print(f"synthetic q0, tol={mp.nstr(tol,2)}:", rel, mp.nstr(res,3) if res else None)
# with cubo q (50-digit value): print to 55 digits
import json
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
S = [mpf(s_) for s_ in J["Sdk"]]
q = S[1]/S[0]
print("cubo q to 55:", mp.nstr(q,55))
print("q - 1/2 - 2ln2:", mp.nstr(q - mpf(1)/2 - 2*mp.log(2), 20))
rel, res = pslq([q, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-45'))
print("cubo q pslq tol 1e-45:", rel, mp.nstr(res,3) if res else None)
