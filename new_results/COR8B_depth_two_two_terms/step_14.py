
import os, signal, time
try:
    os.kill(36962, 9)
    print("killed 36962")
except Exception as e:
    print("kill err:", e)
time.sleep(1)
# verify gone
import subprocess
r = subprocess.run(['/usr/bin/pgrep','-l','-f','step_09.py'],capture_output=True,text=True,timeout=20)
print("remaining:", r.stdout or "(none)")
