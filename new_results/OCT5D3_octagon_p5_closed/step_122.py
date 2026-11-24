import subprocess, sys, os, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
for f in os.listdir("."):
    if f.startswith("gtab2_d360"): os.remove(f)
log = open("verify_run.log","w")
proc = subprocess.Popen([sys.executable, "-u", "verify_oct5d3.py"], stdout=log, stderr=subprocess.STDOUT)
time.sleep(60)
print("poll at 60s: running" if proc.poll() is None else "done rc=%d" % proc.poll())
print(open("verify_run.log").read())
