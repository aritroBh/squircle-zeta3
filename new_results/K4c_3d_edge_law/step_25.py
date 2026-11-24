
test = r'''"""step_04d.py -- cubo outer-refinement test: outer 1/48 & 1/64, inner 1/64."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
cert = [mpf(s) for s in J["Sdk"]]
for HO in [mpf(1)/48]:
    t0=time.time()
    SC, st = chamber_S(NC, 4, HO, mpf(1)/64, mpf('4.5'), mpf('4.5'), HI2=mpf(1)/64)
    el = time.time()-t0
    print(f"outer HO=1/{int(1/HO)}, inner 1/64 ({el:.0f}s) stats={st}")
    for k in range(5):
        print(f"   [d^{k}] err = {mp.nstr(abs(SC[k]-cert[k]),3)}", flush=True)
'''
with open('step_04d.py','w') as f: f.write(test)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_04d.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
