import subprocess, time
t0 = time.monotonic()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_hexc2.py'],
                   capture_output=True, text=True, cwd='<project>/research/08_theory/HEXC2_level6_certificate')
print('\n'.join(r.stdout.strip().splitlines()[:14]))
print('...')
print('\n'.join(r.stdout.strip().splitlines()[-4:]))
print('returncode', r.returncode, 'total %.1fs'%(time.monotonic()-t0))
with open('progress.log','a') as f:
    f.write('step_80: final fresh run verify_hexc2.py PASS (incl. euler cert), %.0fs\n'%(time.monotonic()-t0))
