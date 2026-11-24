import subprocess, sys
r = subprocess.run([sys.executable, "step_hp.py", "w5", "650"], capture_output=True, text=True, timeout=2700)
print(r.stdout); print("ERR", r.stderr[-300:])
r = subprocess.run([sys.executable, "step_hp.py", "m4", "890"], capture_output=True, text=True, timeout=2700)
print(r.stdout); print("ERR", r.stderr[-300:])
