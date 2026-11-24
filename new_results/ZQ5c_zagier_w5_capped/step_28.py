import os
cwd = os.getcwd()
kscript = r'''import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J, cands = w5pipe.load_cands('out/d29', 29)
names=J['names']; d=len(names); N=len(cands)
print(f"d29: N={N} d={d} names={names}", flush=True)
K, comps, sym_dicts, rank = w5pipe.kernel_and_comps(cands, d, 'd29')
json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank, saturated=False),
          open('out/d29/kernel5.json','w'))
print(f"saved out/d29/kernel5.json: kernel {len(K)}, rank {rank}", flush=True)
with open('progress.log','a') as f:
    f.write(f"step_05_d29_kernel: N=77 kernel dim {len(K)}, rank {rank}\n")
'''
open('step_05_d29_kernel.py','w').write(kscript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_05_d29_kernel.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2400)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
