
import os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
t0=time.time()
r = subprocess.run([PY, 'step_10_J192.py'], cwd=FIX, capture_output=True, text=True, timeout=1800)
print(f"({time.time()-t0:.1f}s)")
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else "")
