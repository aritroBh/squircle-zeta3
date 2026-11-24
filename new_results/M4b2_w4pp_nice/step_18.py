import subprocess
base = '<project>'
here = base + '/research/08_theory/M4b2_w4pp_nice'
PY = base + '/.venv/bin/python'
r = subprocess.run(['nice', PY, 'step06_stage2_exact.py'], cwd=here, capture_output=True, text=True)
open(here + '/step06_stage2_exact.out', 'w').write(r.stdout + r.stderr)
print('rc=', r.returncode)
print(r.stdout[-2000:])
if r.returncode: print('STDERR:', r.stderr[-1500:])
