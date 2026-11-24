
import time
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_07.py MODEL (a) FINAL h96: k0 err 4.1e-62 (2sqrt5/3); errs est: k1 1.9e-49, k2 9.7e-48, k3 2.0e-46, k4 1.6e-45 (45+ digits). VALUES: k1=4.310068730207151345980382494723116602408017857946782637185, k2=2.175841230811692460507624858609765605772129518494031707693, k3=-2.330156241024368053549557284389126408582825209655730929953, k4=0.004482674056558090243652254900043497274132120333173616794591. NUMERICAL 45-49 digits.\n")
step08 = r'''"""step_08.py -- MODEL (b): 3 axes + 4 body diagonals, production grids."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S, outer_splits
s3 = 1/mp.sqrt(3)
NB = [[1,0,0],[0,1,0],[0,0,1],[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
print("outer splits:", [mp.nstr(x,8) for x in outer_splits(NB)], flush=True)
res = {}
for tag, HO, HI, HI2 in [("h32", mpf(1)/32, mpf(1)/32, mpf(1)/64), ("h48", mpf(1)/48, mpf(1)/64, mpf(1)/64)]:
    t0=time.time()
    S, st = chamber_S(NB, 4, HO, HI, mpf('4.5'), mpf('4.5'), HI2=HI2)
    el = time.time()-t0
    res[tag] = [mp.nstr(v,70) for v in S]
    print(f"{tag} outer=1/{int(1/HO)} inner=1/{int(1/HI)}/1/{int(1/HI2)} ({el:.0f}s) stats={st}", flush=True)
    for k in range(5):
        print(f"   [d^{k}]S = {mp.nstr(S[k],58)}", flush=True)
    json.dump({"Sdk":res[tag],"stats":st,"HO":str(HO),"HI":str(HI),"HI2":str(HI2)}, open(f"modelB_{tag}.json","w"), indent=1)
# exact S(1) from sympy: -36 + 21sqrt3 + (9/13)sqrt(35515-20504sqrt3)
import sympy as sp
S1 = sp.N(-36 + 21*sp.sqrt(3) + sp.Rational(9,13)*sp.sqrt(35515-20504*sp.sqrt(3)), 60)
print("exact S(1) sympy:", S1)
S48 = [mpf(s) for s in res["h48"]]; S32 = [mpf(s) for s in res["h32"]]
print("k=0 err vs exact:", mp.nstr(abs(S48[0]-mpf(str(S1))),3))
for k in range(5):
    print(f"[d^{k}] |h48-h32| = {mp.nstr(abs(S48[k]-S32[k]),3)}")
'''
with open('step_08.py','w') as f: f.write(step08)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_08.py'], capture_output=True, text=True, timeout=5400)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
