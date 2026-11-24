import subprocess, time
t0=time.time()
r = subprocess.run(["<project>/.venv/bin/python", "verify_srf3r.py"], capture_output=True, text=True, timeout=1500)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
print("wall:", time.time()-t0)
