import subprocess, sys
r = subprocess.run([sys.executable, "step_hp.py", "w4pp", "850"], capture_output=True, text=True, timeout=2700)
print(r.stdout); print("ERR", r.stderr[-300:])
r = subprocess.run([sys.executable, "step_hp.py", "w5", "690"], capture_output=True, text=True, timeout=2700)
print(r.stdout); print("ERR", r.stderr[-300:])
