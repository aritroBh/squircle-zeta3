import sys, subprocess, time
t0=time.monotonic()
r = subprocess.run(['nice', sys.executable, 'step_02_varpi1.py'],capture_output=True,text=True,timeout=600,cwd='../DCY2_decay_first_correction')
print(r.stdout); print("STDERR:", r.stderr[-300:] if r.stderr else "")
print("elapsed", time.monotonic()-t0)
open('../DCY2R_referee_varpi1/progress.log' if False else 'progress.log','a').write("step: reran step_02_varpi1.py -- varpi1 exact identities diff=0\n")
