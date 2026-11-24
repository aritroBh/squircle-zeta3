import subprocess, time, os
dst = '<project>/research/08_theory/M4i_w5_noL'
PY='<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice','-n','10',PY,'lrel_build.py','420'], capture_output=True, text=True, cwd=dst, timeout=1500)
el = time.monotonic()-t0
print('elapsed %.1fs rc=%d' % (el, r.returncode))
print(r.stdout)
print('STDERR:', r.stderr[-3000:])
