import subprocess, time, os
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_zlad1r.py'],
                   capture_output=True, text=True, timeout=1500)
print('exit:', r.returncode, f'({time.time()-t0:.0f}s wall)')
print(r.stdout)
print('STDERR tail:', r.stderr[-800:])
