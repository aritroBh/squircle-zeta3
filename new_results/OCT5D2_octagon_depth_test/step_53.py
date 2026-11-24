import subprocess, os, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
if os.path.exists("gtab_d90_N320.pkl"): os.remove("gtab_d90_N320.pkl")
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run([PY, "verify_oct5d2.py"], capture_output=True, text=True, timeout=900)
print("returncode:", r.returncode, " wall time: %.1fs" % (time.time()-t0))
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
