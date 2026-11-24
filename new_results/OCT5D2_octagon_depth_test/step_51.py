import subprocess, os, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
# remove the cached high-precision G-table so verify builds its own fresh dps-90 table
for f in ["gtab_d90_N320.pkl"]:
    if os.path.exists(f): os.remove(f)
t0 = time.time()
r = subprocess.run(["python3", "verify_oct5d2.py"], capture_output=True, text=True, timeout=600)
print("returncode:", r.returncode, " wall time: %.1fs" % (time.time()-t0))
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
