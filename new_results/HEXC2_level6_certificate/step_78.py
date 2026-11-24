import subprocess, time
t0 = time.monotonic()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_hexc2.py'],
                   capture_output=True, text=True, cwd='<project>/research/08_theory/HEXC2_level6_certificate')
tail = r.stdout.strip().splitlines()[-8:]
print('\n'.join(tail))
print('returncode', r.returncode, 'total %.1fs'%(time.monotonic()-t0))
with open('progress.log','a') as f:
    f.write('step_75: independent fresh re-run verify_hexc2.py PASS, %.0fs\n'%(time.monotonic()-t0))
import os
print(sorted(os.listdir('.')))
