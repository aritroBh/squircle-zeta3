import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import time
import w5lib as W
import pipe
C, cands = pipe.phase_enum(5, B=8)
pipe.print_cands(5, cands, C)
K = pipe.phase_kernel(5, C, cands)
import json
json.dump(dict(K=K), open(pipe.cand_dir(5)+"/kernel.json","w"))
'''
open('steps/step_27_d5_enum.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_27_d5_enum.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_27_d5_enum: D=5 enum+kernel\n")
