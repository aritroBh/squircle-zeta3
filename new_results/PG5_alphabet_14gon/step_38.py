import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_02_s5_ranks.py','345'],capture_output=True,text=True,timeout=1500)
print("elapsed:", time.time()-t0)
print(r.stdout)
print(r.stderr[-2000:] if r.stderr else '')
