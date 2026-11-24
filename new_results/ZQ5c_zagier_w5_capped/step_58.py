import os, subprocess, time
cwd = os.getcwd()
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_20_d13cap_lll.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2700, env={**os.environ,'FULL':'1'})
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-1500:] if r.stderr else "")
