import subprocess, sys, os, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
for f in os.listdir("."):
    if f.startswith("gtab2_d360"): os.remove(f)
t0 = time.time()
r = subprocess.run([sys.executable, "verify_oct5d3.py"], capture_output=True, text=True, timeout=600)
print("returncode:", r.returncode, " wall time: %.1fs" % (time.time()-t0))
print(r.stdout)
print("STDERR tail:", r.stderr[-2000:] if r.stderr else "")
