import subprocess, time
t0 = time.time()
r = subprocess.run(["<project>/.venv/bin/python", "-u", "verify_dod5d.py"],
                   capture_output=True, text=True, timeout=290)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
print("wall: %.0fs, rc=%d" % (time.time()-t0, r.returncode))
