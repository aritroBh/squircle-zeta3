import subprocess, time
t0 = time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_tower.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR:', r.stderr[-1000:])
print('elapsed: %.0fs' % (time.time()-t0))
