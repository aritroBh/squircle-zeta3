"""PSLQ for c_2 over the weight-4 alternating-MZV basis at two precisions. Usage: python id_c2.py"""
import json
from w4_common import *
lam_names = ["z3*lam", "pi^4", "pi^2*lam^2", "lam^4", "Li4(1/2)"]
out = {}
for dps, tolexp in ((360, 330), (345, 315)):
    set_dps(dps)
    c = json.load(open(os.path.join(RES, "id_coeffs_dps360.json")))
    c2 = mp.mpf(c["c"][2])
    lam = mp.log(2); z3 = mp.zeta(3); pi = mp.pi
    vals = [z3*lam, pi**4, pi**2*lam**2, lam**4, mp.polylog(4, mp.mpf(1)/2)]
    rel = mp.pslq([c2] + vals, maxcoeff=10**6, maxsteps=10**7, tol=mp.mpf(10)**(-tolexp))
    print("dps", mp.mp.dps, "tol 1e-%d" % tolexp, "relation (c2, z3 lam, pi^4, pi^2 lam^2, lam^4, Li4(1/2)):", rel)
    if rel:
        res = sum(r*x for r, x in zip(rel, [c2] + vals)); print("   residual", mp.nstr(res, 5), " dps", mp.mp.dps)
    out[dps] = rel
json.dump({str(k): v for k, v in out.items()}, open(os.path.join(RES, "id_c2_pslq.json"), "w"))
# closed form check, dps 360
set_dps(360); lam = mp.log(2)
cf = -(3*mp.pi**4 + 20*mp.pi**2*lam**2 + 20*lam**4 + 120*mp.zeta(3)*lam)/960
c2 = mp.mpf(json.load(open(os.path.join(RES, "id_coeffs_dps360.json")))["c"][2])
print("c2 - closed form:", mp.nstr(c2 - cf, 5), " dps", mp.mp.dps)
set_dps(300); lam = mp.log(2)
cf = -(3*mp.pi**4 + 20*mp.pi**2*lam**2 + 20*lam**4 + 120*mp.zeta(3)*lam)/960
c2 = mp.mpf(json.load(open(os.path.join(RES, "id_coeffs_dps300.json")))["c"][2])
print("c2(dps300 data) - closed form:", mp.nstr(c2 - cf, 5), " dps", mp.mp.dps)
