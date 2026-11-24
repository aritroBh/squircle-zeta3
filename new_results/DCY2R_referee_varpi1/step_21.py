import sys, subprocess, time
print(sys.executable)
t0=time.monotonic()
r = subprocess.run(['nice', sys.executable, 'dcy2_audit.py'],capture_output=True,text=True,timeout=1500,cwd='../DCY2_decay_first_correction')
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.stderr else "")
print("elapsed", time.monotonic()-t0)
