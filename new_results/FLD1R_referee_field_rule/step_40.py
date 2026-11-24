import subprocess, time
t0=time.time()
r = subprocess.run(['nice','<project>/.venv/bin/python','verify_fld1.py'],
                   capture_output=True, text=True, cwd='.', timeout=1500)
print('exit', r.returncode, f'{time.time()-t0:.0f}s')
print(r.stdout[-1800:])
if r.stderr: print('STDERR:', r.stderr[-400:])
open('progress.log','a').write(f'step06: FLD1 verify_fld1.py rerun in referee sandbox: exit {r.returncode} in {time.time()-t0:.0f}s\n')
