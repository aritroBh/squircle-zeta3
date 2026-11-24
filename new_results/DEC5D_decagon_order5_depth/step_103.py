
import subprocess, time
t0 = time.time()
r = subprocess.run(["<project>/.venv/bin/python", "verify_dec5d.py"],
                   capture_output=True, text=True, cwd=".", timeout=290)
print("elapsed %.1fs, returncode %d" % (time.time()-t0, r.returncode))
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
