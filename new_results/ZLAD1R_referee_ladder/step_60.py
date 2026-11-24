import subprocess, time
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_zlad1r.py'],
                   capture_output=True, text=True, timeout=1200)
print('exit:', r.returncode, f'({time.time()-t0:.0f}s)')
print(r.stdout)
if r.returncode != 0:
    print('STDERR:', r.stderr[-1000:])
open('progress.log','a').write(f'step_57: post-interruption re-run of verify_zlad1r.py: exit {r.returncode} ({time.time()-t0:.0f}s)\n')
