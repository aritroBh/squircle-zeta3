import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
mp.mp.dps=60
print("L(chi17,5) =", mp.nstr(vals["L"], 30))
print("T = L/sqrt17 =", mp.nstr(vals["T"], 30))
'''
open('steps/step_54c_d17_num.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_54c_d17_num.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_54c_d17_num: D=17 numerics\n")
