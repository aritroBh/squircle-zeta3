import subprocess, time
t0=time.time()
r = subprocess.run(["<project>/.venv/bin/python","verify_fn2.py"],
                   capture_output=True, text=True, cwd=".", timeout=600)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "(none)")
print("wall:", time.time()-t0)
