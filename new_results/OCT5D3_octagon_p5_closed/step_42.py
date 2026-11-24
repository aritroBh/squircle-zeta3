import subprocess, sys, os, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
log = open("build_big.log","w")
proc = subprocess.Popen([sys.executable, "gengine2p.py", "1650", "5700", "13"],
                        stdout=log, stderr=subprocess.STDOUT)
print("launched PID", proc.pid)
# also build a dps=1000/N=3500 table for the 1000-digit closure (cache separately, cheap) after big one?
# time estimate: poll after 30s
time.sleep(30)
print("poll: running" if proc.poll() is None else "done")
print(open("build_big.log").read())
