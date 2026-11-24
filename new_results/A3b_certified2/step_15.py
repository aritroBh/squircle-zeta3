
import subprocess, sys, time
t0 = time.time()
r = subprocess.run([sys.executable, 'w4_certified.py', '400', '500'],
                   capture_output=True, text=True, timeout=3600)
print("returncode:", r.returncode, " wall: %.1fs" % (time.time()-t0))
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
