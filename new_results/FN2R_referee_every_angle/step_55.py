import subprocess, time
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","verify_fn2r.py"],
                   capture_output=True, text=True, cwd=".", timeout=1500)
out = r.stdout
print(out)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
open("verify_fn2r.out","w").write(out)
print("wall:", time.time()-t0)
