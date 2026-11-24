import subprocess, os
base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
if os.path.exists(here + '/verifier_nice.log'): os.remove(here + '/verifier_nice.log')
PY = base + '/.venv/bin/python'
r = subprocess.run(['nice', PY, 'verify_W4pp_nice.py'], cwd=here, capture_output=True, text=True)
open(here + '/verify_W4pp_nice.out', 'w').write(r.stdout + r.stderr)
print('rc=', r.returncode)
print(r.stdout[-3500:])
if r.returncode: print('STDERR:', r.stderr[-2000:])
