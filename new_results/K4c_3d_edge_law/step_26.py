
import time
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_04c/04d: inner refinement useless (error identical 1/64 vs 1/96) -> outer-limited. OUTER 1/48 + inner 1/64: cubo errs 4.9e-60/3.0e-60/1.4e-58/1.0e-56/4.7e-55 vs K4b cert (k=0..4). Driver PRODUCTION-READY (k4c_driver.py).\n")
step05 = r'''"""step_05.py -- MODEL (a): 6 face-diagonal normals, [d^k]S k=0..4, two grids."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S, outer_splits
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
print("outer splits:", [mp.nstr(x,8) for x in outer_splits(NA)], flush=True)
res = {}
for tag, HO, HI in [("h32", mpf(1)/32, mpf(1)/32), ("h48", mpf(1)/48, mpf(1)/48)]:
    t0=time.time()
    S, st = chamber_S(NA, 4, HO, HI, mpf('4.5'), mpf('4.5'), HI2=HI)
    el = time.time()-t0
    res[tag] = [mp.nstr(v,70) for v in S]
    print(f"{tag} ({el:.0f}s) stats={st}", flush=True)
    for k in range(5):
        print(f"   [d^{k}]S = {mp.nstr(S[k],55)}", flush=True)
    json.dump({"Sdk":res[tag],"stats":st,"HO":str(HO),"HI":str(HI)}, open(f"modelA_{tag}.json","w"), indent=1)
S1exact = 2*mp.sqrt(5)/3
S48 = [mpf(s) for s in res["h48"]]; S32 = [mpf(s) for s in res["h32"]]
print("k=0 err vs 2sqrt5/3:", mp.nstr(abs(S48[0]-S1exact),3))
for k in range(5):
    print(f"[d^{k}] |h48-h32| = {mp.nstr(abs(S48[k]-S32[k]),3)}")
'''
with open('step_05.py','w') as f: f.write(step05)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_05.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
