import subprocess, time
dst = '<project>/research/08_theory/M4i_w5_noL'
PY='<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice','-n','10',PY,'nol_build.py','60'], capture_output=True, text=True, cwd=dst, timeout=900)
print('elapsed %.1fs rc=%d' % (time.monotonic()-t0, r.returncode))
print(r.stdout)
print('STDERR:', r.stderr[-3000:])
