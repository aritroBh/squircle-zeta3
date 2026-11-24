import subprocess, os, sys, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
# kill any stray verify processes
subprocess.run(["pkill", "-f", "verify_oct5d3.py"], capture_output=True)
time.sleep(2)
for f in os.listdir("."):
    if f.startswith("gtab2_d360"): os.remove(f)
t0 = time.time()
r = subprocess.run([sys.executable, "-u", "verify_oct5d3.py"], capture_output=True, text=True, timeout=420)
print("rc:", r.returncode, " wall: %.0fs" % (time.time()-t0))
out = r.stdout
print(out[-3500:])
print("STDERR tail:", r.stderr[-1000:] if r.stderr else "")
