import subprocess, time
t0 = time.time()
r = subprocess.run(["python", "verify_cor3.py"], capture_output=True, text=True, timeout=590)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "(none)")
print(f"return code {r.returncode}, wall time {time.time()-t0:.1f}s")
