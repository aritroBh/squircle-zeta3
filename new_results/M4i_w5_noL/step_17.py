import os, subprocess
print('cwd:', os.getcwd())
PY='<project>/.venv/bin/python'
r = subprocess.run([PY,'-c','import flint, mpmath; print("flint ok", flint.__version__ if hasattr(flint,"__version__") else "?", mpmath.__version__)'], capture_output=True, text=True)
print(r.stdout, r.stderr)
# quick check: how long does nice_search solve take + drop_L as-is
import time
t0=time.monotonic()
r = subprocess.run(['nice','-n','10',PY,'drop_L.py'], capture_output=True, text=True, cwd=os.getcwd())
print('drop_L elapsed %.1fs rc=%d' % (time.monotonic()-t0, r.returncode))
print(r.stdout)
print('STDERR:', r.stderr[-2000:])
