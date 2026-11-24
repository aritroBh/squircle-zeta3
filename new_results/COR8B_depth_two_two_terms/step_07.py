
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_02.py'],capture_output=True,text=True,timeout=110)
print(r.stdout); print(r.stderr); print("elapsed", time.time()-t)
