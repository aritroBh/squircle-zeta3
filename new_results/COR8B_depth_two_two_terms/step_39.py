
import subprocess, os
r = subprocess.run(['/usr/bin/pgrep','-l','-f','COR8B'],capture_output=True,text=True,timeout=20)
print(r.stdout)
mypid = os.getpid()
for line in r.stdout.strip().splitlines():
    pid = int(line.split()[0])
    if pid != mypid:
        try:
            os.kill(pid, 9); print("killed", pid, line[:80])
        except Exception as e:
            print("err", pid, e)
