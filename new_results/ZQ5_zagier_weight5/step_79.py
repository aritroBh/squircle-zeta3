import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import w5lib as W
import pipe
D=21
t0=time.time()
C, cands = pipe.phase_enum(D, B=8)
print(f"enum {time.time()-t0:.0f}s")
pipe.print_cands(D, cands, C)
t0=time.time()
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
print(f"kernel {time.time()-t0:.0f}s dim {len(K)}")
'''
open('steps/step_48_d21_enum.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_48_d21_enum.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_48_d21_enum: D=21 enum+kernel\n")
