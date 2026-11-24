import subprocess, time
t0 = time.time()
r = subprocess.run(['python', 'verify_cor7.py'], capture_output=True, text=True, timeout=600)
print(r.stdout)
print("STDERR tail:", r.stderr[-800:] if r.stderr else "")
print(f"total wall time {time.time()-t0:.0f}s, returncode {r.returncode}")
