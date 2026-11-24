import sys, subprocess, json
J = json.load(open('out/d17/candidates.json'))
print("N candidates:", len(J['cands']))
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
t0=time.time()
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
print(f"kernel {time.time()-t0:.0f}s dim {len(K)}")
'''
open('steps/step_53b_d17_kernel.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_53b_d17_kernel.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_53b_d17_kernel\n")
