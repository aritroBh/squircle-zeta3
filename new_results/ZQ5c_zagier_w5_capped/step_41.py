import os
cwd = os.getcwd()
kscript = r'''import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J = json.load(open('out/d13/candidates_cap.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
print(f"d13-cap: N={N} d={d}", flush=True)
K, comps, sym_dicts, rank = w5pipe.kernel_and_comps(cands, d, 'd13cap')
json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank, saturated=False),
          open('out/d13/kernel5_cap.json','w'))
print(f"saved out/d13/kernel5_cap.json: kernel {len(K)}, rank {rank}", flush=True)
with open('progress.log','a') as f:
    f.write(f"step_14_d13cap_kernel: cap N={N} kernel dim {len(K)}, rank {rank}\n")
'''
open('step_14_d13cap_kernel.py','w').write(kscript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_14_d13cap_kernel.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=1500)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
print("=== BG values log ===")
print(open('step_12_bg.log').read())
