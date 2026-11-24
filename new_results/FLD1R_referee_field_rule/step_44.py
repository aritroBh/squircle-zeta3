import subprocess, time
t0=time.time()
r = subprocess.run(['nice','<project>/.venv/bin/python','verify_fld1r.py'],
                   capture_output=True, text=True, cwd='.', timeout=600)
print('exit', r.returncode, f'{time.time()-t0:.0f}s')
print(r.stdout)
if r.stderr: print('STDERR:', r.stderr[-500:])
open('progress.log','a').write(f'step08: verify_fld1r.py: exit {r.returncode} in {time.time()-t0:.0f}s\n')
