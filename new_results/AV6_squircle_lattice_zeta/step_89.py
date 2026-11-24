code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time, json
setdps(65)
taus = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4]
vals = {}
t0=time.time()
for sgn in [1, -1]:
    for tau0 in taus:
        tau = sgn*tau0
        p_ = 1/(1-tau)
        for s_ in [5, 6]:
            ZC_ = ZC(p_, s_, K=500, R=16)
            vals[(tau, s_)] = str(ZC_)
        print(f"tau={tau:+.2f} done ({time.time()-t0:.0f}s)")
json.dump({f"{a},{b}":v for (a,b),v in vals.items()}, open("taugrid56.json","w"))
'''
exec(code)
with open("progress.log","a") as f: f.write("step_49: tau grid s=5,6\n")
