import subprocess, time
t0=time.time()
r = subprocess.run(['nice','<project>/.venv/bin/python','verify_fld1r.py'],
                   capture_output=True, text=True, cwd='.', timeout=600)
print('exit', r.returncode, f'{time.time()-t0:.0f}s')
print(r.stdout[-900:])
open('progress.log','a').write(f'step11: verify_fld1r.py final reproducibility run: exit {r.returncode} in {time.time()-t0:.0f}s\n')
print('='*72)
print(open('progress.log').read())
