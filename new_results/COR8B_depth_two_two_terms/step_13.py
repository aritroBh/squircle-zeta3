
import subprocess, os
mypid = os.getpid()
print("mypid", mypid)
try:
    r = subprocess.run(['/usr/bin/pgrep','-l','-f','step_'],capture_output=True,text=True,timeout=20)
    print(r.returncode); print(r.stdout); print(r.stderr)
except Exception as e:
    print("pgrep err", e)
