
import os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
t0 = time.time()
r = subprocess.run([PY, 'planner_sweep_ratio.py', '1.0002'], cwd=FIX, capture_output=True, text=True, timeout=3600)
el = time.time()-t0
print(f"ELAPSED {el:.1f}s")
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else "")
with open(os.path.join(WORK,'progress.log'),'a') as f:
    tail = r.stdout.strip().split(chr(10))
    f.write("D13s turn9: sweep ratio 1.0002 (%.0fs): %s\n" % (el, " | ".join(tail[:3])))
