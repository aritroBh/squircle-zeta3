import subprocess, time
t0=time.time()
r = subprocess.run(['<project>/.venv/bin/python','verify_pg6.py'],
                   capture_output=True, text=True, cwd='.', timeout=600)
print(r.stdout)
print("STDERR (tail):", r.stderr[-1500:] if r.stderr else "")
print("elapsed: %.0fs"%(time.time()-t0))
