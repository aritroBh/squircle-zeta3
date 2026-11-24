
step09 = r'''"""step_09.py -- MODEL (b) outer h=1/96 certification run."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s3 = 1/mp.sqrt(3)
NB = [[1,0,0],[0,1,0],[0,0,1],[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
t0=time.time()
S, st = chamber_S(NB, 4, mpf(1)/96, mpf(1)/64, mpf('4.5'), mpf('4.5'), HI2=mpf(1)/64)
el = time.time()-t0
print(f"h96outer ({el:.0f}s) stats={st}", flush=True)
vals = []
for k in range(5):
    v = S[k]
    if mp.im(v) != 0:
        assert abs(mp.im(v)) < mpf('1e-100')
        v = mp.re(v)
    vals.append(mp.nstr(v,70))
    print(f"   [d^{k}]S = {mp.nstr(v,58)}", flush=True)
json.dump({"Sdk":vals,"stats":st,"HO":"1/96","HI":"1/64","HI2":"1/64"}, open("modelB_h96.json","w"), indent=1)
J48 = json.load(open("modelB_h48.json")); S48=[mpf(s) for s in J48["Sdk"]]
J32 = json.load(open("modelB_h32.json")); S32=[mpf(s) for s in J32["Sdk"]]
S96 = [mpf(s) for s in vals]
print("trend:")
for k in range(5):
    d1 = abs(S48[k]-S32[k]); d2 = abs(S96[k]-S48[k])
    r = d1/d2 if d2>0 else mpf('inf')
    est = d2/(mp.sqrt(r)-1) if r>1 else d2   # conservative: one halving-step gain sqrt(r)
    print(f"[d^{k}] |h48-h32|={mp.nstr(d1,3)} |h96-h48|={mp.nstr(d2,3)} ratio={mp.nstr(r,3)} -> h96 err est~{mp.nstr(est,3)}")
'''
with open('step_09.py','w') as f: f.write(step09)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_09.py'], capture_output=True, text=True, timeout=5400)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
