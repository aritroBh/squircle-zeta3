
import time
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_05.py MODEL (a) [6 face diagonals = tetrakis hexahedron]: k=0 = 1.490711984999859797606115... matches 2sqrt5/3 to 4.0e-62 (PROVES S(1)=2sqrt5/3; step_02 float value was hull-noise). h48: k1=4.310068730207151..., k2=2.175841230811692..., k3=-2.330156241024368..., k4=0.004482674056558090...; h48-h32 diffs: k1 5.8e-31, k2 1.7e-29, k3 2.0e-28. NUMERICAL pending h64 trend.\n")
step06 = r'''"""step_06.py -- MODEL (a) h=1/64 for convergence trend."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
t0=time.time()
S, st = chamber_S(NA, 4, mpf(1)/64, mpf(1)/64, mpf('4.5'), mpf('4.5'), HI2=mpf(1)/64)
el = time.time()-t0
print(f"h64 ({el:.0f}s) stats={st}", flush=True)
for k in range(5):
    print(f"   [d^{k}]S = {mp.nstr(S[k],55)}", flush=True)
json.dump({"Sdk":[mp.nstr(v,70) for v in S],"stats":st}, open("modelA_h64.json","w"), indent=1)
J48 = json.load(open("modelA_h48.json")); S48=[mpf(s) for s in J48["Sdk"]]
J32 = json.load(open("modelA_h32.json")); S32=[mpf(s) for s in J32["Sdk"]]
print("convergence trend (diffs):")
for k in range(5):
    d1 = abs(S48[k]-S32[k]); d2 = abs(S[k]-S48[k])
    r = d1/d2 if d2>0 else mpf('inf')
    print(f"[d^{k}] |h48-h32|={mp.nstr(d1,3)}  |h64-h48|={mp.nstr(d2,3)}  ratio={mp.nstr(r,3)}")
'''
with open('step_06.py','w') as f: f.write(step06)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_06.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
