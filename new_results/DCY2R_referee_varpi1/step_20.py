import subprocess, time
t0=time.monotonic()
r = subprocess.run(['bash','-c','cd ../DCY2_decay_first_correction && nice <home>/Downloads/"Zeta(3)"/.venv/bin/python dcy2_audit.py 2>&1 | head -60'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print(r.stderr[-500:] if r.stderr else "")
print("elapsed", time.monotonic()-t0)
open('progress.log','a').write("step: reran dcy2_audit.py (N=400,2000,10000) -- see output\n")
