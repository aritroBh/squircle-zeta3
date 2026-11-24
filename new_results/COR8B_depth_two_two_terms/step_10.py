
import subprocess, signal, os, time
# try pkill via /bin/kill? first try pgrep alternatives; sandbox may block ps. Try pkill directly.
for pat in ['step_02.py', 'step_01.py']:
    try:
        r = subprocess.run(['/usr/bin/pkill','-9','-f',pat],capture_output=True,text=True,timeout=20)
        print(pat, r.returncode, r.stdout, r.stderr)
    except Exception as e:
        print(pat, 'ERR', e)
time.sleep(1)
# quick self-test of speed again
import mpmath as mp
mp.mp.dps = 50
t=time.time()
print(mp.polylog(2, mp.mpf('0.5')), time.time()-t)
