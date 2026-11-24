import subprocess, os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
PY = base + '/.venv/bin/python'
import time
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'verify_town1.py'], capture_output=True, text=True, cwd=workdir)
print("stdout:")
print(r.stdout)
print("stderr:", r.stderr[-2000:] if r.stderr else "")
print("rc:", r.returncode, f"({time.time()-t0:.1f}s)")
