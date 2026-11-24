import subprocess, time, sys
t0 = time.time()
r = subprocess.run([sys.executable, 'verify_cor7.py'], capture_output=True, text=True, timeout=900)
print(r.stdout)
print("STDERR tail:", r.stderr[-800:] if r.stderr else "")
print(f"total wall time {time.time()-t0:.0f}s, returncode {r.returncode}")
