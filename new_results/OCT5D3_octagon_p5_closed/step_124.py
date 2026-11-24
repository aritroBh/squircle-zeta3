import subprocess, sys, os, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
for f in os.listdir("."):
    if f.startswith("gtab2_d360"): os.remove(f)
log = open("verify_run.log","w")
proc = subprocess.Popen([sys.executable, "-u", "verify_oct5d3.py"], stdout=log, stderr=subprocess.STDOUT)
t0=time.time()
while proc.poll() is None and time.time()-t0 < 420:
    time.sleep(15)
rc = proc.poll()
print("rc:", rc, " wall: %.0fs" % (time.time()-t0))
print(open("verify_run.log").read())
if rc is None:
    proc.kill(); print("KILLED - too slow")
