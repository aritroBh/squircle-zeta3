import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json
import w5lib as W
import pipe
C, cands = pipe.phase_enum(13, B=8)
pipe.print_cands(13, cands, C)
K = pipe.phase_kernel(13, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(13)+"/kernel.json","w"))
print("kernel dim", len(K))
'''
open('steps/step_43_d13_enum.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_43_d13_enum.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_43_d13_enum: D=13 enum+kernel\n")
