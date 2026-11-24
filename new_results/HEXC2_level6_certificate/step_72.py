import subprocess, time
t0 = time.monotonic()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_hexc2.py'],
                   capture_output=True, text=True, cwd='<project>/research/08_theory/HEXC2_level6_certificate')
print(r.stdout)
print('STDERR:', r.stderr[-2000:] if r.stderr else '')
print('returncode', r.returncode, 'total %.1fs'%(time.monotonic()-t0))
