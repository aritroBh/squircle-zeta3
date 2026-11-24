import os, time, subprocess
pid=6401
try:
    os.kill(pid,0); print('process alive')
except OSError as e:
    print('process check:', e)
r=subprocess.run(['bash','-lc','ls /proc/6401/ 2>/dev/null | head; cat /proc/6401/stat 2>/dev/null | cut -d" " -f3,14,15'],capture_output=True,text=True)
print(r.stdout, r.stderr)
