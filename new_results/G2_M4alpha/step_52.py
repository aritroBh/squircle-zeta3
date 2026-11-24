import subprocess, sys, time
cwd = "<project>/research/08_theory/G2_M4alpha"
t0 = time.time()
r = subprocess.run([sys.executable, "step_07.py"], capture_output=True, text=True, cwd=cwd, timeout=1200)
print("time", time.time() - t0)
print(r.stdout)
print("STDERR", r.stderr[-500:])
