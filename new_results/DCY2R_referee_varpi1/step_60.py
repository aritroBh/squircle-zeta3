import subprocess, sys, time
t0=time.monotonic()
r = subprocess.run(['nice', sys.executable, 'verify_dcy2r.py'],capture_output=True,text=True,timeout=900)
print(r.stdout)
print("STDERR:", r.stderr[-500:] if r.stderr else "")
print("total wall:", time.monotonic()-t0)
