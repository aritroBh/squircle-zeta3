import subprocess, sys, time
cwd = "<project>/research/08_theory/G2_M4alpha"
t0 = time.time()
r = subprocess.run([sys.executable, "step_05.py", "640"], capture_output=True, text=True, cwd=cwd, timeout=1800)
print("time", time.time() - t0)
print(r.stdout[-1000:])
print("STDERR", r.stderr[-300:])
